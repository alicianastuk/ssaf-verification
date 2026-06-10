"""SSAF verification suite runner. Executes every test_*.py; exits nonzero on any failure."""
import glob, subprocess, sys

failures = []
for f in sorted(glob.glob("test_*.py")):
    print(f"=== {f} ===")
    r = subprocess.run([sys.executable, f])
    if r.returncode != 0:
        failures.append(f)
print("\n" + ("SUITE PASS" if not failures else f"SUITE FAIL: {failures}"))
sys.exit(1 if failures else 0)
