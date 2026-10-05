#!/usr/bin/env python3
"""Live read-only TUI for a PRRS run: counters, position on the PES, bond events.

Design constraints, deliberately:
- Pure observer. Nothing inside runs/ is ever written; the full-resolution PES
  snapshot goes to .scratch/prrs_watch/<run>/pes_live.png. A monitor that could
  perturb a run would be worse than no monitor.
- Fail-closed display. A value that is not recorded shows as "n/a (reason)",
  never as zero; a map that cannot be drawn states why instead of drawing
  something that looks like a result.
- No new dependencies. numpy is already a project dependency; matplotlib is
  already used by scripts/p3_node_bundle/p3_figure.py. Without matplotlib the
  text panels still work and the map says why it is absent.
- The PES panel reuses the p3_figure.py style constants (viridis, white
  contours, meV-above-minimum colour scale) so the live view and the paper
  figure cannot drift apart visually.

Data sources, all produced by the running search anyway:
- runs/<run>/network.json  -- rewritten continuously by publish(); status,
  outcome_counts, budgets, attempts, reactions, ts_candidates. May be absent
  for probe runs; that is reported, not invented.
- runs/<run>/trials/*/observations.jsonl -- appended per frame with a flush
  (src/prrs/io.py append_json). Parsed incrementally from a byte offset, so a
  poll costs only the new lines.
- trials/*/result.json -- completion marker; its absence is what "RUN" means.

Usage:
  python3 scripts/prrs_watch.py --run runs/p2_round6
  python3 scripts/prrs_watch.py --run runs/p3_x --grid runs/p3_pes2d [--barrier 316]
  python3 scripts/prrs_watch.py --demo                 # synthetic data, no run
  python3 scripts/prrs_watch.py --run ... --once       # one text frame + PNG, no TUI

Keys: q quit, +/- halve/double the text refresh interval.
"""

import argparse
import collections
import glob
import json
import os
import shutil
import sys
import termios
import time
from collections import namedtuple

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import patheffects as pe
    HAVE_MPL = True
except ImportError:  # text panels keep working; the map reports the gap
    HAVE_MPL = False

Frame = namedtuple("Frame", "phase step time_fs epot fmax unc tracked events")

# palette for text, RGB tuples (truecolor) or quantised to 256 by the emitter
WHITE = (222, 228, 234)
DIM = (138, 148, 158)
GREEN = (126, 222, 126)
RED = (255, 108, 108)
YELLOW = (255, 210, 96)
CYAN = (110, 202, 224)
BG_DARK = (16, 20, 24)
RAMP = " .:-=+*#%@"


def load_json(path):
    """utf-8-sig because the search itself loads network.json with that preamble."""
    try:
        with open(path, "r", encoding="utf-8-sig") as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def fmt(value, digits=3, unit=""):
    return "n/a" if value is None else f"{value:.{digits}f}{unit}"


def hhmmss(seconds):
    seconds = max(0, int(seconds))
    return f"{seconds // 3600}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


def bar(frac, width=10):
    frac = 0.0 if frac < 0 else min(1.0, frac)
    filled = int(round(frac * width))
    return "█" * filled + "░" * (width - filled)


# ---------------------------------------------------------------- data layer


class Trial:
    """One trials/<id>/ directory; observations.jsonl parsed incrementally."""

    def __init__(self, path, tid, sink):
        self.path = path
        self.tid = tid
        self.sink = sink  # run-level deque for bond events
        self.obs_path = os.path.join(path, "observations.jsonl")
        self.result_path = os.path.join(path, "result.json")
        self.offset = 0
        self.partial = ""
        self.frames = []
        self.last = None
        self.done = False
        self.classification = None
        self.mtime = 0.0
        # Set once result.json exists -- search writes it only after the trial ends --
        # and observations.jsonl has been read to its end; final_stamp is that file's
        # (inode, size) at the time, None if the trial left no observations.
        self.finished = False
        self.final_stamp = None

    @property
    def n_frames(self):
        return len(self.frames)

    def poll(self):
        try:
            stat = os.stat(self.obs_path)
        except OSError:
            stat = None
        stamp = None if stat is None else (stat.st_ino, stat.st_size)
        if self.finished:
            if stamp == self.final_stamp:
                return  # finished and untouched: one stat, no reads
            # the directory was replaced (a resume rolls trials back and reuses ids). A
            # bare stat can see that up to a minute late over NFS; for a rare event on a
            # display that is the right price, an open per finished trial per poll is not.
            self.finished, self.done, self.classification = False, False, None
            if stamp is None or self.final_stamp is None or stamp[0] != self.final_stamp[0]:
                self.offset, self.partial, self.frames, self.last = 0, "", [], None
        result = load_json(self.result_path)
        if result is not None:
            self.done = True
            self.classification = result.get("classification")
        if stat is None:
            if self.done:
                self.finished, self.final_stamp = True, None
            return
        self.mtime = stat.st_mtime
        if stat.st_size < self.offset:  # replaced or truncated: start over
            self.offset, self.partial, self.frames, self.last = 0, "", [], None
        try:
            with open(self.obs_path, "rb") as handle:
                handle.seek(self.offset)
                raw = handle.read()
        except OSError:
            return
        if not raw:
            if self.done:
                self.finished, self.final_stamp = True, stamp
            return
        self.offset += len(raw)
        text = self.partial + raw.decode("utf-8", errors="replace")
        lines = text.split("\n")
        self.partial = lines.pop()  # trailing half-line stays buffered
        for line in lines:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except ValueError:
                continue  # a torn line is skipped, never half-interpreted
            events = record.get("events") or {}
            frame_no = record.get("frame")
            for kind in ("formed", "broken"):
                for pair in events.get(kind) or []:
                    self.sink.append((time.time(), self.tid, frame_no, kind, pair))
            self.frames.append(Frame(
                record.get("phase"), record.get("step"), record.get("time_fs"),
                record.get("potential_eV"), record.get("max_force_eV_A"),
                record.get("uncertainty_eV_A"), record.get("tracked") or {},
                events))
        if self.frames:
            self.last = self.frames[-1]


class RunData:
    """Polls network.json + trials/; keeps last good values on transient NFS hiccups."""

    def __init__(self, run_dir):
        self.run_dir = run_dir
        self.name = os.path.basename(os.path.normpath(run_dir))
        self.network = None
        self.status = "no network.json yet"
        self.trials = {}
        self.events = collections.deque(maxlen=200)
        self.started = None  # approx wall clock of the oldest observation
        self.started_done = False
        self.network_stamp = None

    def poll(self):
        # publish() replaces network.json atomically, so a new version always has a new
        # inode; re-parse only then. The stamp moves only on a good parse, so a torn or
        # unreadable read is retried on the next poll. open + fstat rather than stat:
        # the run is usually written on another host, and over NFS only an open
        # revalidates, while a bare stat can serve cached attributes for up to a minute.
        path = os.path.join(self.run_dir, "network.json")
        try:
            with open(path, "rb") as handle:
                stat = os.fstat(handle.fileno())
            stamp = (stat.st_ino, stat.st_mtime_ns, stat.st_size)
        except OSError:
            stamp = None
        if stamp is not None and stamp != self.network_stamp:
            net = load_json(path)
            if net is not None:
                self.network = net
                self.status = net.get("status", "?")
                self.network_stamp = stamp
        try:
            names = os.listdir(os.path.join(self.run_dir, "trials"))
        except OSError:
            names = []
        for name in sorted(names):
            path = os.path.join(self.run_dir, "trials", name)
            if not os.path.isdir(path):
                continue
            trial = self.trials.get(name)
            if trial is None:
                trial = Trial(path, name, self.events)
                self.trials[name] = trial
            trial.poll()
            if not self.started_done and trial.mtime:
                self.started = trial.mtime if self.started is None \
                    else min(self.started, trial.mtime)
        self.started_done = True

    def active(self):
        running = [t for t in self.trials.values() if not t.done and t.last]
        return max(running, key=lambda t: t.mtime) if running else None

    def counters(self):
        net = self.network or {}
        trials = list(self.trials.values())
        reactions = net.get("reactions")
        ts_candidates = net.get("ts_candidates")
        return {
            "trials_total": len(trials),
            "trials_done": sum(1 for t in trials if t.done),
            "frames_total": sum(t.n_frames for t in trials),
            "attempts": net.get("attempts", "n/a"),
            "reactions": len(reactions) if isinstance(reactions, list) else "n/a",
            "ts_candidates": len(ts_candidates) if isinstance(ts_candidates, list) else "n/a",
            "budgets": net.get("budgets") or {},
            "outcomes": net.get("outcome_counts") or {},
        }


class DemoRun(RunData):
    """Synthetic walk between two wells; exercises the full pipeline with no run."""

    KX, KY = "q_PT", "r_OO"

    def __init__(self):
        super().__init__("demo")
        self.status = "running (DEMO -- synthetic)"
        self.network = {
            "status": "running",
            "attempts": 6,
            "reactions": [{"a": 1}, {"a": 2}],
            "ts_candidates": [1, 2, 3, 4, 5],
            "budgets": {"chemical_nodes_used": 3, "chemical_nodes_available": 8,
                        "microstates_total": 4},
            "outcome_counts": {"completed": 10, "ts_candidate": 2, "rejected": 2},
        }
        self.started = time.time() - 4267
        self.walk_q, self.walk_r = self._make_walk()
        self.pos = 0
        for i in range(14):
            tid = f"s{i:02d}L0300_t000000"
            self.trials[tid] = self._fake_trial(tid, done=True)
        self.trials["s14L0300_t000000"] = self._fake_trial(
            "s14L0300_t000000", done=False)

    @staticmethod
    def _make_walk():
        t = np.linspace(0.0, 1.0, 400)
        rng = np.random.default_rng(7)
        q = -0.62 + 0.60 * t ** 1.6 + 0.006 * rng.standard_normal(len(t))
        r = 2.545 - 0.19 * t ** 1.2 + 0.004 * rng.standard_normal(len(t))
        return q, r

    def _fake_trial(self, tid, done):
        trial = Trial.__new__(Trial)
        trial.tid, trial.sink, trial.done = tid, self.events, done
        trial.classification = ["elastic", "reaction", None][hash(tid) % 3] if done else None
        trial.frames, trial.last = [], None
        trial.mtime = time.time()
        trial.offset, trial.partial = 0, ""
        trial.obs_path, trial.result_path = f"demo:{tid}", f"demo:{tid}"
        rng = np.random.default_rng(hash(tid) % (2 ** 32))
        well = -0.55 if hash(tid) % 2 else 0.60
        for k in range(int(rng.integers(20, 45))):
            q = well + 0.06 * rng.standard_normal()
            r = 2.50 + 0.03 * rng.standard_normal()
            trial.frames.append(Frame("quenched", k, 0.1 * k, None, None, None,
                                      {self.KX: float(q), self.KY: float(r)}, {}))
        if trial.frames:
            trial.last = trial.frames[-1]
        return trial

    def poll(self):
        self.pos = min(self.pos + 2, len(self.walk_q) - 1)
        trial = self.trials["s14L0300_t000000"]
        frame = Frame("free", self.pos, 0.1 * self.pos, None, None, None,
                      {self.KX: float(self.walk_q[self.pos]),
                       self.KY: float(self.walk_r[self.pos])}, {})
        if self.pos == 180:
            self.events.append((time.time(), trial.tid, 180, "broken", [3, 7]))
        trial.frames.append(frame)
        trial.last = frame
        trial.mtime = time.time()


# ------------------------------------------------------------------ PES map


def load_grid(path):
    """pes2d npz or a directory of pes2d_relaxed_shard*.npz, merged as p3_merge does.

    Returns (grid, None) or (None, reason); an empty surface is reported, never
    drawn as flat.
    """
    base = os.path.expanduser(path)
    files = []
    if os.path.isdir(base):
        merged = os.path.join(base, "pes2d_relaxed.npz")
        files = [merged] if os.path.exists(merged) \
            else sorted(glob.glob(os.path.join(base, "pes2d_relaxed_shard*.npz")))
    elif os.path.exists(base):
        files = [base]
    if not files:
        return None, f"no pes2d grid found at {path}"
    Q = R = E = conv = None
    for fname in files:
        blob = np.load(fname)
        if Q is None:
            Q, R = blob["Q"], blob["R"]
            E = np.full((len(R), len(Q)), np.nan)
            conv = np.full((len(R), len(Q)), -1, dtype=np.int8)
        if "rows" in blob:  # shard: owns strided full rows, scatter them in
            for i in blob["rows"]:
                E[i] = blob["E"][i]
                conv[i] = blob["converged"][i] if "converged" in blob else 1
        else:  # already merged: the whole grid is present
            E[:] = blob["E"]
            conv[:] = blob["converged"] if "converged" in blob else 1
    E = np.where(conv == 0, np.nan, E)  # unconverged points are not surface values
    if not np.isfinite(E).any():
        return None, "grid present but 0 points finite -- treating as absent"
    return {"Q": Q, "R": R, "E": (E - np.nanmin(E)) * 1e3}, None


def area_to_cells(buf, w_cells, h_cells):
    """RGBA buffer -> (2h, w) RGB area average -> half-block cells (top, bottom)."""
    h, w, _ = buf.shape
    yi = np.round(np.linspace(0, h, 2 * h_cells + 1)).astype(int)
    xi = np.round(np.linspace(0, w, w_cells + 1)).astype(int)
    small = np.empty((2 * h_cells, w_cells, 3))
    for j in range(2 * h_cells):
        for i in range(w_cells):
            block = buf[yi[j]:yi[j + 1], xi[i]:xi[i + 1]]
            small[j, i] = block.reshape(-1, 3).mean(axis=0)
    return [[(tuple(int(v) for v in small[2 * r, c]),
              tuple(int(v) for v in small[2 * r + 1, c]))
             for c in range(w_cells)] for r in range(h_cells)]


class PesMap:
    def __init__(self, grid_arg, cap_mev, barrier_mev):
        self.grid, self.grid_error = (None, None)
        if grid_arg:
            self.grid, self.grid_error = load_grid(grid_arg)
        self.cap = cap_mev
        self.barrier = barrier_mev

    def render(self, data, active, kx, ky, w_cells, h_cells, png_path):
        """Returns (cells, title_note, problem). cells: rows of (top_rgb, bot_rgb)."""
        if not HAVE_MPL:
            return None, None, "matplotlib missing -- map unavailable (text panels fine)"
        if self.grid_error:
            problem = self.grid_error
        if w_cells < 12 or h_cells < 7:
            return None, None, "panel too small for the map"
        ghost_x, ghost_y, ghost_e = [], [], []
        for trial in data.trials.values():
            if trial is active:
                continue  # the active trial is the trail, not ghosts
            for fr in trial.frames:
                x, y = fr.tracked.get(kx), fr.tracked.get(ky)
                if x is not None and y is not None:
                    ghost_x.append(x)
                    ghost_y.append(y)
                    ghost_e.append(fr.epot)
        trail = []
        if active is not None:
            for fr in active.frames:
                x, y = fr.tracked.get(kx), fr.tracked.get(ky)
                if x is not None and y is not None:
                    trail.append((x, y, fr.phase))
        problem = None if self.grid is not None else self.grid_error
        if self.grid is None and len(ghost_x) + len(trail) < 2:
            problem = problem or "no frames with tracked coordinates yet"

        dark = [c / 255 for c in BG_DARK]
        fig = plt.figure(figsize=(w_cells / 10.0, 2.0 * h_cells / 10.0), dpi=100)
        fig.patch.set_facecolor(dark)
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        ax.set_facecolor(dark)
        for spine in ax.spines.values():
            spine.set_color("#55606a")
            spine.set_linewidth(0.8)
        ax.set_xticks([])
        ax.set_yticks([])

        if self.grid is not None:
            g = self.grid
            levels = np.linspace(0.0, self.cap, 19)
            ax.contourf(g["Q"], g["R"], np.clip(g["E"], 0, self.cap), levels=levels,
                        cmap="viridis", extend="max")
            ax.contour(g["Q"], g["R"], g["E"], levels=np.arange(50, self.cap, 50),
                       colors="white", linewidths=0.45, alpha=0.55)
            if self.barrier and 0 < self.barrier < self.cap:
                ax.contour(g["Q"], g["R"], g["E"], levels=[self.barrier],
                           colors="white", linewidths=1.3, alpha=0.9)
            ax.set_xlim(g["Q"][0], g["Q"][-1])
            ax.set_ylim(g["R"][0], g["R"][-1])
            note = f"relaxed grid, 0..{self.cap:.0f} meV (cap)"
        else:
            all_x = np.array(ghost_x + [p[0] for p in trail], dtype=float)
            all_y = np.array(ghost_y + [p[1] for p in trail], dtype=float)
            energies = [e * 1e3 for e in ghost_e if e is not None]
            if len(all_x) >= 2:
                if energies:
                    ax.scatter(all_x[:len(ghost_x)], all_y[:len(ghost_y)], c=energies,
                               cmap="viridis", s=9, lw=0, alpha=0.85,
                               vmin=min(energies), vmax=max(energies) + 1e-6)
                else:
                    ax.scatter(all_x[:len(ghost_x)], all_y[:len(ghost_y)], s=9,
                               lw=0, alpha=0.6, color="#9aa4ae")
            ax.set_xlim(*self._padded(all_x))
            ax.set_ylim(*self._padded(all_y))
            note = f"frames only, meV above run min (n={len(ghost_x) + len(trail)})"

        if trail:
            tx = np.array([p[0] for p in trail])
            ty = np.array([p[1] for p in trail])
            quench = np.array([p[2] in ("quench", "quenched") for p in trail])
            halo = [pe.Stroke(linewidth=3.4, foreground="white"), pe.Normal()]
            free = ~quench
            if free.any():
                ax.plot(tx[free], ty[free], "-", color="#ff4444", lw=2.0, alpha=0.9,
                        solid_capstyle="round", path_effects=halo, zorder=6)
            if quench.any():
                ax.plot(tx[quench], ty[quench], ".", color="#60aaff", ms=4, zorder=6)
            ax.plot(tx[-1], ty[-1], "o", ms=10, mfc="#ff4444", mec="white", mew=1.6,
                    zorder=8)
            ax.axhline(ty[-1], color="#ffd24d", lw=0.6, alpha=0.5, zorder=5)
            ax.axvline(tx[-1], color="#ffd24d", lw=0.6, alpha=0.5, zorder=5)

        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
        fig.savefig(png_path, dpi=200, facecolor=fig.get_facecolor())
        plt.close(fig)
        return area_to_cells(buf, w_cells, h_cells), note, problem

    @staticmethod
    def _padded(values):
        if len(values) == 0 or not np.isfinite(values).any():
            return 0.0, 1.0
        lo, hi = np.nanmin(values), np.nanmax(values)
        if hi - lo < 1e-9:
            hi = lo + 1.0
        pad = 0.05 * (hi - lo)
        return lo - pad, hi + pad


# ------------------------------------------------------------------ terminal


def detect_colors(mode, env):
    if mode != "auto":
        return mode
    colorterm = env.get("COLORTERM", "")
    if "truecolor" in colorterm or "24bit" in colorterm:
        return "truecolor"
    term = env.get("TERM", "")
    if "256color" in term or "kitty" in term or "wezterm" in term:
        return "256"
    return "mono"


def rgb256(rgb):
    r, g, b = rgb
    return 16 + 36 * int(round(r / 255 * 5)) + 6 * int(round(g / 255 * 5)) \
        + int(round(b / 255 * 5))


class Screen:
    """Hand-rolled ANSI screen: no curses, so truecolor half-block cells are trivial."""

    def __init__(self, colors, once=False):
        self.colors = colors
        self.once = once
        self.cols, self.rows = self._size()
        self.cells = [[[" ", None, None, False] for _ in range(self.cols)]
                      for _ in range(self.rows)]
        self._saved_attrs = None

    @staticmethod
    def _size():
        cols, rows = shutil.get_terminal_size((100, 32))
        return max(40, cols), max(14, rows)

    def resize(self):
        cols, rows = self._size()
        if (cols, rows) != (self.cols, self.rows):
            self.cols, self.rows = cols, rows
            self.cells = [[[" ", None, None, False] for _ in range(cols)]
                          for _ in range(rows)]
            self._prev = None

    def enter(self):
        if self.once:
            return
        attrs = termios.tcgetattr(0)
        self._saved_attrs = attrs
        raw = termios.tcgetattr(0)
        raw[3] &= ~(termios.ECHO | termios.ICANON)
        raw[6][termios.VMIN] = 0
        raw[6][termios.VTIME] = 0
        termios.tcsetattr(0, termios.TCSANOW, raw)
        sys.stdout.write("\x1b[?1049h\x1b[?25l\x1b[2J")

    def leave(self):
        if self.once:
            return
        if self._saved_attrs is not None:
            termios.tcsetattr(0, termios.TCSANOW, self._saved_attrs)
        sys.stdout.write("\x1b[0m\x1b[?25h\x1b[?1049l")
        sys.stdout.flush()

    def clear(self):
        for row in self.cells:
            for cell in row:
                cell[0], cell[1], cell[2], cell[3] = " ", None, None, False

    def put(self, x, y, text, fg=None, bold=False):
        if y < 0 or y >= self.rows:
            return
        for i, ch in enumerate(str(text)):
            if x + i < 0 or x + i >= self.cols:
                continue
            self.cells[y][x + i] = [ch, fg, None, bold]

    def blit(self, x, y, cells):
        for j, row in enumerate(cells):
            if y + j >= self.rows:
                break
            for i, (top, bot) in enumerate(row):
                if x + i >= self.cols:
                    break
                if self.colors == "mono":
                    lum = int((np.mean(top) + np.mean(bot)) / 2)
                    ch = RAMP[min(len(RAMP) - 1, lum * len(RAMP) // 256)]
                    self.cells[y + j][x + i] = [ch, DIM, None, False]
                else:
                    self.cells[y + j][x + i] = ["▀", tuple(top), tuple(bot), False]

    def box(self, x, y, w, h, title=None, fg=DIM):
        if w < 4 or h < 2 or x + w > self.cols or y + h > self.rows:
            return
        self.put(x, y, "┌" + "─" * (w - 2) + "┐", fg)
        for yy in range(y + 1, y + h - 1):
            self.put(x, yy, "│", fg)
            self.put(x + w - 1, yy, "│", fg)
        self.put(x, y + h - 1, "└" + "─" * (w - 2) + "┘", fg)
        if title:
            self.put(x + 2, y, f" {title} "[: w - 4], fg)

    def render(self):
        """Differential redraw: only changed cell runs are emitted -- a full
        truecolor frame over SSH costs hundreds of KB, a diff costs a few."""
        if self.once:
            for row in self.cells:
                print("".join(cell[0] for cell in row).rstrip())
            return
        prev = getattr(self, "_prev", None)
        if prev is None or len(prev) != self.rows or len(prev[0]) != self.cols:
            prev = [[None] * self.cols for _ in range(self.rows)]
        out = []
        for r in range(self.rows):
            row, prow = self.cells[r], prev[r]
            c = 0
            while c < self.cols:
                if row[c] == prow[c]:
                    c += 1
                    continue
                out.append(f"\x1b[{r + 1};{c + 1}H")
                style = None  # terminal style state is unknown after a cursor move
                while c < self.cols and row[c] != prow[c]:
                    ch, fg, bg, bold = row[c]
                    state = (fg, bg, bold)
                    if state != style:
                        seq = "\x1b[0m"
                        if bold:
                            seq += "\x1b[1m"
                        if self.colors == "truecolor":
                            if fg:
                                seq += f"\x1b[38;2;{fg[0]};{fg[1]};{fg[2]}m"
                            if bg:
                                seq += f"\x1b[48;2;{bg[0]};{bg[1]};{bg[2]}m"
                        elif self.colors == "256":
                            if fg:
                                seq += f"\x1b[38;5;{rgb256(fg)}m"
                            if bg:
                                seq += f"\x1b[48;5;{rgb256(bg)}m"
                        out.append(seq)
                        style = state
                    out.append(ch)
                    c += 1
            # deep-copy the cells: compose() mutates the live buffer in place,
            # so sharing inner lists would make every later frame look unchanged
            prev[r] = [[cell[0], cell[1], cell[2], cell[3]] for cell in row]
        out.append("\x1b[0m")
        self._prev = prev
        sys.stdout.write("".join(out))
        sys.stdout.flush()


# ------------------------------------------------------------------- layout


def layout(cols, rows):
    """Shared geometry for compose() and the map renderer's pixel budget."""
    trials_h = min(12, max(6, (rows - 6) // 3))
    mid_top = 5
    mid_h = rows - 6 - trials_h
    two_col = cols >= 104 and mid_h >= 12
    if two_col:
        mw = max(40, min(cols - 44, 78))
        return {"trials_h": trials_h, "mid_top": mid_top, "mid_h": mid_h,
                "two_col": True, "mw": mw, "mh": mid_h,
                "rx": mw + 1, "rw": cols - mw - 1}
    mh = max(10, mid_h - 9)
    return {"trials_h": trials_h, "mid_top": mid_top, "mid_h": mid_h,
            "two_col": False, "mw": cols, "mh": mh, "rx": 0, "rw": cols}


def now_panel(data, active):
    if active is None or active.last is None:
        return [(f"no running trial (status: {data.status})", DIM)]
    fr = active.last
    lines = [(f"trial {active.tid}   RUN", YELLOW)]
    lines.append((f"phase {fr.phase}  step {fr.step}  frame {active.n_frames}", WHITE))
    lines.append((f"t {fmt(fr.time_fs, 1, ' fs')}   E_pot {fmt(fr.epot)} eV", WHITE))
    lines.append((f"Fmax {fmt(fr.fmax)} eV/A   unc {fmt(fr.unc)}", WHITE))
    for key, value in list(fr.tracked.items())[:4]:
        shown = f"{value:.4f}" if isinstance(value, (int, float)) else str(value)
        lines.append((f"{key} = {shown}", CYAN))
    return lines


def event_lines(data, limit):
    out = []
    for ts, tid, frame, kind, pair in sorted(data.events, reverse=True)[:limit]:
        stamp = time.strftime("%H:%M:%S", time.localtime(ts))
        if isinstance(pair, (list, tuple)) and len(pair) == 2:
            bond = f"{pair[0]}-{pair[1]}"
        else:
            bond = str(pair)
        out.append((f"{stamp}  f{frame}  {bond} {kind}  {tid}",
                    RED if kind == "broken" else GREEN))
    return out or [("none observed yet", DIM)]


def trial_lines(data, limit):
    rows = []
    for trial in sorted(data.trials.values(), key=lambda t: t.mtime, reverse=True)[:limit]:
        last = trial.last
        status = "done" if trial.done else "RUN"
        phase = last.phase if last else "-"
        epot = f"{last.epot:.2f}" if last and last.epot is not None else "n/a"
        fg = YELLOW if not trial.done else (WHITE if trial.classification else DIM)
        rows.append((f"{trial.tid}  {status}  {trial.classification or '-'}  "
                     f"{trial.n_frames:>4}f  {phase}  E {epot}", fg))
    return rows or [("no trials yet", DIM)]


def compose(scr, data, map_cells, map_note, map_problem, kx, ky, args, colors_mode,
            png_path):
    cols, rows = scr.cols, scr.rows
    geo = layout(cols, rows)
    scr.clear()
    demo = " [DEMO]" if isinstance(data, DemoRun) else ""
    scr.put(0, 0, f"PRRS WATCH  run={data.name}{demo}", WHITE, bold=True)
    wall = hhmmss(time.time() - data.started) if data.started else "n/a"
    scr.put(0, 1, f"status={data.status}  ~wall {wall}  refresh {args.interval:.1f}s"
                  f"  colors {colors_mode}", DIM)
    cnt = data.counters()
    scr.put(0, 2, f"trials {cnt['trials_done']}/{cnt['trials_total']}  "
                  f"frames {cnt['frames_total']}  attempts {cnt['attempts']}  "
                  f"reactions {cnt['reactions']}  ts_candidates {cnt['ts_candidates']}",
             WHITE)
    budgets = cnt["budgets"]
    avail = budgets.get("chemical_nodes_available")
    if avail:
        used = budgets.get("chemical_nodes_used", 0)
        outcomes = "  ".join(f"{k} {v}" for k, v in sorted(cnt["outcomes"].items()))
        scr.put(0, 3, f"nodes {used}/{avail} {bar(used / max(1, avail))}  "
                      f"microstates {budgets.get('microstates_total', 'n/a')}  {outcomes}",
                 CYAN)
    else:
        scr.put(0, 3, "no network.json (probe run) -- counters from trials only", DIM)

    active = data.active()
    map_title = f"PES  {kx} x {ky}"
    if map_note:
        map_title += f"  --  {map_note}"
    scr.box(0, geo["mid_top"], geo["mw"], geo["mh"], map_title)
    if map_cells is not None:
        scr.blit(1, geo["mid_top"] + 1, map_cells)
    else:
        reason = map_problem or "rendering..."
        scr.put(2, geo["mid_top"] + max(1, geo["mh"] // 2), f"[ {reason} ]", YELLOW)

    now = now_panel(data, active)
    trials_top = rows - geo["trials_h"] - 1
    if geo["two_col"]:
        nh = min(len(now) + 2, geo["mid_h"] - 3)
        scr.box(geo["rx"], geo["mid_top"], geo["rw"], nh, "NOW")
        for i, (text, fg) in enumerate(now):
            scr.put(geo["rx"] + 1, geo["mid_top"] + 1 + i, text[: geo["rw"] - 2], fg)
        ev_top = geo["mid_top"] + nh
        ev_h = geo["mid_top"] + geo["mid_h"] - ev_top
        scr.box(geo["rx"], ev_top, geo["rw"], ev_h, "bond events (newest first)")
        evs = event_lines(data, max(1, ev_h - 2))
        for i, (text, fg) in enumerate(evs):
            scr.put(geo["rx"] + 1, ev_top + 1 + i, text[: geo["rw"] - 2], fg)
    else:
        ny = geo["mid_top"] + geo["mh"]
        half = (cols - 2) // 2
        scr.box(0, ny, half + 1, 9, "NOW")
        for i, (text, fg) in enumerate(now[:7]):
            scr.put(1, ny + 1 + i, text[: half - 1], fg)
        scr.box(half + 1, ny, cols - half - 1, 9, "bond events (newest first)")
        for i, (text, fg) in enumerate(event_lines(data, 7)):
            scr.put(half + 2, ny + 1 + i, text[: cols - half - 3], fg)

    scr.box(0, trials_top, cols, geo["trials_h"],
            f"trials (newest first, {cnt['trials_total']})")
    for i, (text, fg) in enumerate(trial_lines(data, geo["trials_h"] - 2)):
        scr.put(1, trials_top + 1 + i, text[: cols - 4], fg)
    scr.put(0, rows - 1, f"png: {png_path}   q quit  +/- interval", DIM)


# --------------------------------------------------------------------- main


def resolve_axes(data, args):
    if args.x and args.y:
        return args.x, args.y, None
    for trial in data.trials.values():
        for fr in trial.frames:
            keys = [k for k, v in (fr.tracked or {}).items() if v is not None]
            if len(keys) >= 2:
                return (args.x or keys[0]), (args.y or keys[1]), None
    return None, None, "no frame with two tracked coordinates yet (pass --x/--y)"


def main(argv=None):
    ap = argparse.ArgumentParser(description="read-only live watcher for a PRRS run")
    ap.add_argument("--run", help="runs/<run> directory")
    ap.add_argument("--grid", help="pes2d npz or shard directory as the map background")
    ap.add_argument("--x", help="tracked key for the map x axis (default: first tracked)")
    ap.add_argument("--y", help="tracked key for the map y axis (default: second tracked)")
    ap.add_argument("--interval", type=float, default=1.0, help="text refresh seconds")
    ap.add_argument("--map-interval", type=float, default=2.0, help="map refresh seconds")
    ap.add_argument("--cap", type=float, default=450.0, help="grid colour cap in meV")
    ap.add_argument("--barrier", type=float, help="meV: one heavier contour at this level")
    ap.add_argument("--colors", choices=["auto", "truecolor", "256", "mono"], default="auto")
    ap.add_argument("--out-root", default=os.path.join(".scratch", "prrs_watch"))
    ap.add_argument("--png-name", default="pes_live.png")
    ap.add_argument("--once", action="store_true", help="one text frame + PNG, no TUI")
    ap.add_argument("--demo", action="store_true", help="synthetic data, no run needed")
    args = ap.parse_args(argv)

    if not args.demo and not args.run:
        ap.error("--run is required (or use --demo)")
    data = DemoRun() if args.demo else RunData(args.run)
    data.poll()

    out_dir = os.path.join(args.out_root, data.name)
    os.makedirs(out_dir, exist_ok=True)
    png_path = os.path.join(out_dir, args.png_name)

    kx, ky, why_not = resolve_axes(data, args)
    map_obj = PesMap(args.grid, args.cap, args.barrier)
    colors_mode = detect_colors(args.colors, os.environ)
    scr = Screen(colors_mode, once=args.once)

    if args.once:
        geo = layout(scr.cols, scr.rows)
        if kx and ky:
            cells, note, problem = map_obj.render(
                data, data.active(), kx, ky, geo["mw"] - 2, geo["mh"] - 2, png_path)
        else:
            cells, note, problem = None, None, why_not
        compose(scr, data, cells, note, problem, kx or "?", ky or "?", args,
                colors_mode, png_path)
        scr.render()
        print(f"\nmap png: {png_path}" + (f"  ({problem})" if problem else ""))
        return 0

    scr.enter()
    try:
        last_map, map_cells, map_note, map_problem = 0.0, None, None, None
        while True:
            data.poll()
            if kx is None or ky is None:
                kx, ky, why_not = resolve_axes(data, args)
            now = time.time()
            if kx and ky and now - last_map >= args.map_interval:
                geo = layout(scr.cols, scr.rows)
                TRAIL_OWNER_ACTIVE = data.active()
                map_cells, map_note, map_problem = map_obj.render(
                    data, TRAIL_OWNER_ACTIVE, kx, ky, geo["mw"] - 2, geo["mh"] - 2,
                    png_path)
                last_map = now
            scr.resize()
            compose(scr, data, map_cells, map_note, map_problem, kx or "?", ky or "?",
                    args, colors_mode, png_path)
            scr.render()
            t0 = time.time()
            while time.time() - t0 < args.interval:
                try:
                    chunk = sys.stdin.buffer.read1(16)
                except (OSError, ValueError):
                    chunk = b""
                if chunk and chunk[:1] in (b"q", b"Q", b"\x03"):
                    return 0
                if chunk[:1] in (b"+", b"="):
                    args.interval = max(0.2, args.interval / 2)
                elif chunk[:1] == b"-":
                    args.interval = min(5.0, args.interval * 2)
                time.sleep(0.05)
    except KeyboardInterrupt:
        return 0
    finally:
        scr.leave()


if __name__ == "__main__":
    sys.exit(main())
