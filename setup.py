#!/usr/bin/env python
#
# This file is part of g_mmpbsa
#
# Author: Rajendra Kumar
#
#
# g_mmpbsa is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# g_mmpbsa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with g_mmpbsa.  If not, see <http://www.gnu.org/licenses/>.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#============================================================================

# Always prefer setuptools over distutils
from setuptools import setup, Extension, find_packages
from setuptools.extension import Library
from setuptools.command.build_ext import build_ext, customize_compiler
from setuptools.command.install import install
from distutils.command.build import build
import sys
import setuptools
import os
import glob

sys.path.append(os.path.dirname(__file__))
here = os.path.abspath(os.path.dirname(__file__))

gromacs_flags = None
apbs_flags = None
extensions = None 
        
def check_gromacs_dirs():
    ''' Check for GROMACS directories and flags
    '''
    import pkgconfig
    
    # If gromacs not in pkgconfig, return from here
    if not pkgconfig.exists('libgromacs'):
        return
    
    out = dict()
    
    cppflags_t = pkgconfig.re.split('\s+', pkgconfig.cflags('libgromacs'))
    out['cppflags'] = []
    out['ldflags'] = []
    out['include'] = []
    out['lib_dirs'] = []
    out['libs'] = []
    
    # Extract include directory and CXXFLAGS
    for flags in cppflags_t:
        if '-I' in flags:
            out['include'].append(flags[2:])
        else:
            out['cppflags'].append(flags)
            
    # Extract lib directory and LDFLAGS
    ldflags_t = pkgconfig.re.split('\s+', pkgconfig.libs('libgromacs'))
    for flags in ldflags_t:
        if '-L' in flags:
            out['lib_dirs'].append(flags[2:])
        elif '-l' in flags:
            out['libs'].append(flags[2:])
        else:
            out['ldflags'].append(flags)
 
    if '-fopenmp' not in out['ldflags']:
        out['ldflags'].append('-fopenmp')

    
    return out

def include_gromacs_source_headers():
    global  gromacs_flags
    if 'GMX_SRC' not in os.environ:
        raise LookupError('GMX_SRC environment variable not found...')
        
    gmx_src = os.path.join(os.environ['GMX_SRC'], 'src')
    gromacs_flags['include'].append(gmx_src)
    for dir_name in glob.glob(f'{gmx_src}/gromacs/*/include'):
        gromacs_flags['include'].append(dir_name)
    

def extract_gromacs_flags():
    ''' Extract gromacs include, lib and other flags for compilation
    '''
    global  gromacs_flags
    
    # At first check if gromacs is already available in standard path
    if 'GMX_INSTALL' not in os.environ:
        gromacs_flags = check_gromacs_dirs()
        if gromacs_flags is not None:
            raise LookupError('Gromacs package not found... Use GMX_INSTALL environment variable to provide GROMACS path.')
    else:
        # If gromacs is not available at standard path check for GMX_INSTALL environment variable,
        # add it to pkg-config and then extract GROMACS directories and flags
        gmx_install = os.environ['GMX_INSTALL']
        if not os.path.isdir(gmx_install):
            raise LookupError('GROMACS directory {0} not exist...'.format(gmx_install))
        # Check lib name: it could be lib or lib64
        lib_dir = None
        for entry in os.listdir(gmx_install):
            if 'lib' in entry:
                lib_dir = entry
                break
        os.environ['PKG_CONFIG_PATH'] = os.path.join(gmx_install, lib_dir, 'pkgconfig')
        gromacs_flags = check_gromacs_dirs()
    if gromacs_flags is None:
        raise LookupError("gromacs package not found")
        
    include_gromacs_source_headers()

def populate_apbs_flags():
    global apbs_flags
    if not 'APBS_INSTALL' in os.environ:
        raise LookupError('APBS_INSTALL environment variable not found...')
    apbs_install = os.environ['APBS_INSTALL']
    if not os.path.isdir(apbs_install):
        raise LookupError('APBS directory {0} not exist...'.format(apbs_install))
    
    # check lib64 or lib dir exit
    lib_dir = None
    for entry in os.listdir(apbs_install):
        if 'lib' in entry:
            lib_dir = entry
            break

    apbs_flags = dict()
    apbs_flags['include'] = [os.path.join(apbs_install, 'include')]
    apbs_flags['lib_dirs'] = [os.path.join(apbs_install, lib_dir)]
    apbs_flags['ldflags'] = [
        '-lapbs_mg',
        '-lapbs_fem',
        '-lapbs_pmgc',
        '-lapbs_generic',
        '-lapbs_routines',
        '-lmc',
        '-lpunc',
        '-lmaloc',
        '-lvf2c',
        '-lcgcode',
        '-lsuperlu',
    ]

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

def get_extensions():
    '''Extensions that are need to be compiled
    '''
    global extensions
    extensions = [ Extension(
        'g_mmpbsa.g_mmpbsa',
        [   'src/pywrapper.cpp',
            'src/apbs.cpp',
            'src/mmpbsa.cpp',
            'src/energy2bfac.cpp',
            'src/nsc.cpp'
            ],
        include_dirs=[ get_pybind_include(), get_pybind_include(user=True), 
                      'src', ] + gromacs_flags['include'] + apbs_flags['include'],
        extra_compile_args = ['-fopenmp'],
        library_dirs=gromacs_flags['lib_dirs'] + apbs_flags['lib_dirs'],
        libraries=gromacs_flags['libs'],
        runtime_library_dirs = gromacs_flags['lib_dirs']  + apbs_flags['lib_dirs'],
        language='c++',
        extra_link_args= gromacs_flags['ldflags'] + apbs_flags['ldflags'],
        ),]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++17'):
        return '-std=c++17'
    elif has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'linux':
        c_opts['unix'] += [ '-static-libstdc++'] # Got From https://github.com/pypa/manylinux/issues/118

    def build_extensions(self):
        # Check for -stdlib=libc++ on macos-clang
        if sys.platform == 'darwin':
            # first check ""-stdlib=libstdc++" is available,
            # if available means gcc is used in place of clang
            if has_flag(self.compiler, '-stdlib=libstdc++'):
                self.c_opts['unix'] += ['-stdlib=libstdc++']
                
            # Only in case of clang, so check for this flag
            elif has_flag(self.compiler, '-stdlib=libc++'):
                self.c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args += opts
            if sys.platform == 'linux':
                ext.extra_link_args += ['-static-libstdc++']

        # Remove "-Wstrict-prototypes" flags
        customize_compiler(self.compiler)
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass

        build_ext.build_extensions(self)

populate_apbs_flags()
extract_gromacs_flags()
get_extensions()
setup(
    ext_modules=extensions,
    cmdclass={'build_ext': BuildExt},
    packages=find_packages(),
    include_package_data=True,
)
