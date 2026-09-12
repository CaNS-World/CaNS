#!/usr/bin/env python3
"""Run the built-in Poisson correctness check on valid boundary/grid layouts."""
import itertools
import os
from pathlib import Path
import re
import signal
import subprocess

# Fixture, grid, process grid, pencil axis, distributed tridiagonal solve.
LAYOUTS = [
    ('even', (8, 10, 12), (1, 1), 1, False),
    ('odd', (9, 11, 13), (1, 2), 1, False),
    ('odd', (9, 11, 13), (2, 2), 2, False),
    ('odd', (9, 11, 13), (2, 3), 1, True),
    ('odd', (9, 11, 13), (3, 2), 2, True),
    ('odd', (9, 11, 13), (1, 2), 3, False),
    ('stretched', (9, 11, 9), (1, 1), 1, False),
    ('stretched', (9, 11, 9), (1, 2), 1, True),
    ('thin', (2, 5, 7), (2, 1), 2, False),
    ('thin', (5, 2, 7), (2, 1), 1, False),
    ('thin', (5, 7, 2), (1, 2), 1, False),
]
BCS = [(pair,)*3 for pair in ('PP', 'DD', 'NN', 'DN', 'ND')]
BCS += [('DD', 'NN', 'DN'), ('ND', 'DN', 'PP')]
INVERSE = dict(PP='PP', DD='NN', NN='DD', DN='ND', ND='DN')


def main():
    testdir = Path(__file__).resolve().parent
    rundir = testdir.parent.parent/'run'
    exe = rundir/'cans'
    output = rundir/'poisson'
    os.environ.setdefault('OMPI_MCA_rmaps_base_oversubscribe', '1')
    os.environ.setdefault('OMP_NUM_THREADS', '2')
    for index, (layout, bc) in enumerate(itertools.product(LAYOUTS, BCS), 1):
        fixture, grid, dims, axis, dtdma = layout
        name = f'{index:02d}-{fixture}-{"-".join(bc)}-{dims[0]}x{dims[1]}-p{axis}-d{int(dtdma)}'
        case = output/name
        (case/'data').mkdir(parents=True, exist_ok=True)
        text = (testdir/f'input-{fixture}.nml').read_text()
        updates = {
            'ng(1:3)': ', '.join(map(str, grid)),
            'dims(1:2)': f'{dims[0]}, {dims[1]}, ipencil_axis = {axis}',
            'is_poisson_dtdma': 'T' if dtdma else 'F',
            'cbcpre(0:1,1:3)': ', '.join(repr(c) for pair in bc for c in pair),
        }
        for component in range(1, 4):
            updates[f'cbcvel(0:1,1:3,{component})'] = ', '.join(repr(c) for pair in bc for c in INVERSE[pair])
        for key, value in updates.items():
            text = re.sub(r'^'+re.escape(key)+r'\s*=.*$', key+' = '+value, text, flags=re.M)
        (case/'input.nml').write_text(text)
        with (case/'run.log').open('w') as log:
            process = subprocess.Popen(['mpirun', '-n', str(dims[0]*dims[1]), str(exe)], cwd=case,
                                       stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                status = process.wait(timeout=60)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                raise RuntimeError(f'MPI timeout: {case}/run.log') from None
        log = (case/'run.log').read_text()
        if status != 0 or '*** Fim ***' not in log or 'ERROR:' in log:
            raise RuntimeError(f'Failed {case}\n{log[-3000:]}')
        print(f'PASS {name}', flush=True)
    print(f'PASS: {len(LAYOUTS)*len(BCS)} Poisson cases; logs in {output}')


if __name__ == '__main__':
    main()
