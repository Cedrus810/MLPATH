
import json, time, torch
x = torch.zeros(8, device="cuda")
torch.cuda.synchronize()
start = time.perf_counter()
for _ in range(2000):
    x = x + 1.0
torch.cuda.synchronize()
print(json.dumps({"calls": 2000, "seconds": time.perf_counter() - start}))
