/*
 * This file is part of g_mmpbsa
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2018-2025  Rajendra Kumar
 *
 * g_mmpbsa is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * g_mmpbsa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with g_mmpbsa.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

 #include <pybind11/pybind11.h>
 #include <pybind11/stl.h>
 
 #include <string>
 #include <vector>
 #include <cstdlib>
 #include <cstdio>
 #include <unistd.h>
 #include <signal.h>
 
 
 #include "gromacs/commandline/cmdlineinit.h"
 #include "gromacs/utility/baseversion.h"
 
 namespace py = pybind11;
 
 int apbs(int argc,char *argv[]);
 int mmpbsa(int argc,char *argv[]);
 int energy2bfac (int argc,char *argv[]);
 
 void exit_handler(int s){
     printf(" Caught KeyboardInterrupt in C++\n");
     exit(1);
 }
 
 void register_ctrl_c_signal() {
     struct sigaction sigIntHandler;
     sigIntHandler.sa_handler = exit_handler;
     sigemptyset(&sigIntHandler.sa_mask);
     sigIntHandler.sa_flags = 0;
     sigaction(SIGINT, &sigIntHandler, NULL);
 }
 
 template<typename F>
 void wrapped_gmx_function(std::vector<std::string> argument_vector, F *func) {
     /* Acquire GIL before calling Python code */
     py::gil_scoped_acquire acquire;
 
     char *argv[argument_vector.size()];
     for(size_t n =0; n<argument_vector.size(); n++)
         argv[n] = &argument_vector.at(n)[0];
     
     func(argument_vector.size(), argv);
 }
 
 void wrap_gmx_programs(py::module &m) {
     register_ctrl_c_signal(); // register Ctrl+C for keyboard interruption
     
     std::function<void(std::vector<std::string>)> wrapped_apbs = [](std::vector<std::string>  argument_vector) { 
         wrapped_gmx_function(argument_vector, &apbs);
     };
     
     std::function<void(std::vector<std::string>)> wrapped_mmpbsa = [](std::vector<std::string>  argument_vector) { 
         wrapped_gmx_function(argument_vector, &mmpbsa);
     };
     
     std::function<void(std::vector<std::string>)> wrapped_energy2bfac = [](std::vector<std::string>  argument_vector) { 
         wrapped_gmx_function(argument_vector, &energy2bfac);
     };
     
     m.def("gmx_version", &gmx_version);
     m.def("apbs", wrapped_apbs, py::call_guard<py::gil_scoped_release>());
     m.def("mmpbsa", wrapped_mmpbsa, py::call_guard<py::gil_scoped_release>());
     m.def("energy2bfac", wrapped_energy2bfac, py::call_guard<py::gil_scoped_release>());
 }
 
 
 PYBIND11_MODULE(g_mmpbsa, m) 
 {
    wrap_gmx_programs(m);
 }
 