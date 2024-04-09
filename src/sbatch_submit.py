#!/usr/bin/env python

"""
Submit a batch job to SLURM.
"""

import os
import subprocess
import tempfile
import argparse

def main():
    # Read in the template file
    lib_dir = os.path.dirname(os.path.realpath(__file__))
    with open(lib_dir + '/_sbatch_template.sh') as f:
        template = f.readlines()
    template = ''.join(template)

    # Make the arg parser
    desc = 'Submit code as a dftb job'
    parser = argparse.ArgumentParser(description=desc)

    # Arg specs: Short name, long name, description, default, other args
    args = [['i', 'input', 'Shell cmd to run', None, {'required': True}],
            ['s', 'setup', 'Shell cmd that sets up the environment', ':', {}],
            ['n', 'ntasks', 'Number of tasks', '1', {}],
            ['o', 'nodes', 'Number of nodes', '1', {}],
            ['c', 'cpus', 'Number of CPUs per task', '1', {}],
            ['t', 'time', 'Maximum time to let the job run', '5:0:0', {}],
            ['m', 'mem', 'Memory allocation', '10G', {}],
            ['j', 'jobname', 'jobname', 'test', {}],
            ['f', 'file', 'save to sh file', 'test.sh', {}],
            ['d', 'dir', 'Directory to save logfiles', '.', {}],
            ['q', 'qos', 'Select a QOS on BlueBear', 'bbdefault',
             {'choices': ['bbdefault', 'bbshort']}]]
    argnames = [a[1] for a in args]
    for arg in args:
        parser.add_argument('-' + arg[0],
                            '--' + arg[1],
                            help=arg[2],
                            default=arg[3],
                            **arg[4])
    args = parser.parse_args()

    # Replace the Slurm parameters in the sbatch script
    for argname in argnames:
        template = template.replace('_' + argname, getattr(args, argname))
        if(argname == 'file'):
            sbatchfile=getattr(args, argname)
            #print(sbatchfile)
    # Write a temporary sbatch script
    tf = tempfile.NamedTemporaryFile(mode='w')
    tf.writelines(template)
    tf.flush()

    # Run the job
    #p = subprocess.call('sleep 2 && sbatch ' + tf.name, shell=True)
    p = subprocess.call('cat ' + tf.name +' >run_'+sbatchfile+'.sh', shell=True)


if __name__ == "__main__":
    main()
