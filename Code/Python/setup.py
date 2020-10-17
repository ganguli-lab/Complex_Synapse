# -*- coding: utf-8 -*-
"""
Created on Fri May  4 18:22:44 2018

@author: Subhy
"""
import os.path as osp

try:
    import setuptools
    assert setuptools
except ImportError:
    pass
from distutils.sysconfig import get_python_inc

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from numpy.distutils.misc_util import get_info as get_misc_info
from numpy.distutils.misc_util import get_numpy_include_dirs
from numpy.distutils.system_info import get_info as get_sys_info

# =========================================================================


def in_src(name):
    """Add src/ to file name"""
    return osp.join('src', name)


# =========================================================================
config = Configuration(package_name='complex_synapse')

inc_dirs = [get_python_inc()]
if inc_dirs[0] != get_python_inc(plat_specific=1):
    inc_dirs.append(get_python_inc(plat_specific=1))
inc_dirs.append(get_numpy_include_dirs())
inc_dirs.append('src')

lapack_info = get_sys_info('lapack_opt', 0)  # and {}
npymath_info = get_misc_info("npymath")
all_info = {k: lapack_info[k] + npymath_info[k] for k in lapack_info.keys()}
rearrange = in_src('rearrange_data.c')
MODULE_LOC = 'identify.'

# =============================================================================
config.add_extension(MODULE_LOC + '_baum_welch',
                     sources=[in_src('baum_welch.c'), rearrange],
                     include_dirs=inc_dirs,
                     extra_info=all_info)
# =============================================================================
if __name__ == '__main__':
    setup(**config.todict())
# =============================================================================
