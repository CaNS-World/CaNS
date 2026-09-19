#!/usr/bin/env python3
"""Spatial and temporal convergence of the 2D Taylor-Green vortex (double precision)."""
import os
from pathlib import Path
import re
import signal
import subprocess

import numpy as np


def norms(error):
    return np.array([(np.sqrt(np.mean(e*e)), np.max(abs(e))) for e in error])


def run(n, steps, axis):
    testdir = Path(__file__).resolve().parent
    rundir = testdir.parent.parent/'run'
    grid = [n, n, n]
    grid[axis] = 1
    lengths = [2*np.pi]*3
    lengths[axis] /= n
    dims = ((1, 2), (1, 1), (2, 1))[axis]
    case = rundir/'tgv'/f'{"xyz"[axis]}-n{n}-steps{steps}'
    (case/'data').mkdir(parents=True, exist_ok=True)
    text = (testdir/'input.nml').read_text()
    for key, value in {
        'ng(1:3)': ', '.join(map(str, grid)),
        'l(1:3)': ', '.join(f'{length:.17g}' for length in lengths),
        'dims(1:2)': f'{dims[0]}, {dims[1]}, ipencil_axis = 1',
        'inivel': f"'tgv-2d-{'xyz'[axis]}'",
        'cfl': f'0.5, dtmax = 1., dt_f = {1/steps:.17g}',
        'nstep': f'{steps}, time_max = 1., tw_max = 1.',
        'icheck': f'0, iout0d = 0, iout1d = 0, iout2d = 0, iout3d = 0, isave = {steps}',
    }.items():
        text = re.sub(r'^'+re.escape(key)+r'\s*=.*$', key+' = '+value, text, flags=re.M)
    (case/'input.nml').write_text(text)
    with (case/'run.log').open('w') as log:
        process = subprocess.Popen(['mpirun', '-n', str(dims[0]*dims[1]), str(rundir/'cans')], cwd=case,
                                   stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            status = process.wait(timeout=120)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            raise RuntimeError(f'MPI timeout: {case}/run.log') from None
    log = (case/'run.log').read_text()
    if status != 0 or '*** Fim ***' not in log or 'ERROR:' in log:
        raise RuntimeError(f'Failed {case}\n{log[-3000:]}')
    fields = []
    for field in ('u', 'v', 'w'):
        data = np.fromfile(case/'data'/f'fld_{field}.bin', dtype=np.float64)
        if data.size != np.prod(grid)+2 or not np.isfinite(data).all():
            raise ValueError(f'Invalid field or wrong precision: {case}, {field}')
        time = data[-2]
        if abs(time-1.) > 1e-12 or data[-1] != steps:
            raise ValueError(f'Wrong final time/step: {case}, {field}')
        fields.append(data[:-2].reshape(grid, order='F'))
    if np.max(abs(fields[axis])) > 1e-10:
        raise ValueError(f'Nonzero velocity in homogeneous direction: {case}')
    a, b = (axis+1) % 3, (axis+2) % 3
    return np.array([f.transpose(a, b, axis)[:, :, 0] for f in (fields[a], fields[b])]), time


def check_orders(label, errors, lower, upper):
    errors = np.array(errors)
    if not np.isfinite(errors).all() or not (errors > 0).all():
        raise ValueError(f'Invalid {label} errors: {errors}')
    orders = np.log2(errors[:-1]/errors[1:])
    for order in orders:
        print(f'{label} order: ' + ' '.join(f'{p:.3f}' for p in order.ravel()), flush=True)
    if not ((orders >= lower) & (orders <= upper)).all():
        raise ValueError(f'{label} orders outside [{lower}, {upper}]')


def main():
    os.environ.setdefault('OMP_NUM_THREADS', '2')
    # Allow oversubscription with both Open MPI 4 (ORTE) and 5 (PRRTE).
    os.environ.setdefault('OMPI_MCA_rmaps_base_oversubscribe', '1')
    os.environ.setdefault('PRTE_MCA_rmaps_default_mapping_policy', ':oversubscribe')
    for axis in range(3):
        a, b = 'uvw'[(axis+1) % 3], 'uvw'[(axis+2) % 3]
        print(f'{"XYZ"[axis]} homogeneous:       {a} L2       {a} Linf         {b} L2       {b} Linf', flush=True)
        errors = []
        for n in (16, 32, 64):
            uv, time = run(n, 100, axis)
            xc = (np.arange(n)+.5)*2*np.pi/n
            xf = (np.arange(n)+1.)*2*np.pi/n
            exact = np.array([np.cos(xf[:, None])*np.sin(xc[None, :]),
                              -np.sin(xc[:, None])*np.cos(xf[None, :])])
            error = norms(uv-exact*np.exp(-.2*time))  # nu=0.1, T=1
            errors.append(error)
            print(f'space N={n:2d}:    ' + ' '.join(f'{e:12.5e}' for e in error.ravel()), flush=True)
        check_orders('spatial', errors, 1.9, 2.1)
        # Differences at fixed N cancel the spatial error floor; measure explicit RK3.
        solutions = [run(32, steps, axis)[0] for steps in (20, 40, 80, 160)]
        differences = [norms(a-b) for a, b in zip(solutions, solutions[1:])]
        for steps, error in zip((20, 40, 80), differences):
            print(f'time difference dt={1/steps:g}: ' + ' '.join(f'{e:12.5e}' for e in error.ravel()), flush=True)
        check_orders('temporal', differences, 2.7, 3.3)
    print('PASS: 21 cases; 2D Taylor-Green spatial and temporal convergence in X/Y/Z')


if __name__ == '__main__':
    main()
