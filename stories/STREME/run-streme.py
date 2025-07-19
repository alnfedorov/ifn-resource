import time
from multiprocessing import cpu_count
from subprocess import Popen

import ld

workload = []
for setup in ld.streme.comparisons.glob(f"PLS/*/{ld.streme.target}"):
    targets = setup
    background = setup.parent / ld.streme.background
    saveto = setup.parent / ld.streme.saveto
    workload.append((targets, background, saveto))

running = []
jobs = cpu_count()
while workload:
    while len(running) < jobs and workload:
        targets, background, saveto = workload.pop()
        running.append(Popen([
            "streme",
            "--p", targets, "--n", background,
            "-oc", saveto,
            "--verbosity", "1",
            "--dna",
            "--order", str(ld.streme.markov_order),
            "--seed", str(ld.streme.seed),
            "--minw", str(ld.streme.minw), "--maxw", str(ld.streme.maxw),
            "--hofract", str(ld.streme.hofract),
            "--objfun", "de",
            # "--nmotifs", str(ld.streme.nmotifs),
            "--thresh", str(ld.streme.thresh)
        ]))
    time.sleep(2)

    unfinished = []
    for p in running:
        poll = p.poll()
        if poll is None:
            unfinished.append(p)
        else:
            assert poll == 0, "streme failed"
    running = unfinished
    print(f"{len(workload)} remaining")

for p in running:
    assert p.wait() == 0, "streme failed"
