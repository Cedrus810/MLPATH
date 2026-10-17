# Preliminary probe attempt

The synthetic provenance omitted `electronic_state`. The existing seed export
guard correctly refused before export with:

`ValueError: delivery-gap: no electronic state in provenance; a refinement seed must name the charge/spin of the PES it was found on (re-run pars select)`

This is a probe-construction error, not an algorithm finding or successful
reproduction. Earlier geometry files remain as partial artifacts. The corrected
probe writes a fresh `seed-stages-v2` directory with explicit synthetic charge/spin.
