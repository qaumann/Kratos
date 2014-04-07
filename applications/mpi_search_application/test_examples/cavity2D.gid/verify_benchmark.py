from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
kratos_benchmarking_path = '../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

Text = " "

#
print("verifying  test_fractstep_cont_laplacian.py...")
successful,Msg = benchmarking.RunBenchmark("test_fractstep_cont_laplacian.py", "fractstep_cont_laplacian_benchmarking_ref.txt")

if(successful==True):
    Text += "OK\n"
    print("test_fractstep_cont_laplacian example successful")
else:
    Text += "FAILED\n"
    Text += Msg
    Text += "\n\n"
    print("test_fractstep_cont_laplacian example FAILED")

#
print("verifying  test_fractstep_discrete_laplacian.py...")
successful,Msg = benchmarking.RunBenchmark("test_fractstep_discrete_laplacian.py", "fractstep_discrete_laplacian_benchmarking_ref.txt")


if(successful==True):
    Text += "OK\n"
    print("test_fractstep_discrete_laplacian example successful")
else:
    Text += "FAILED\n"
    Text += Msg
    Text += "\n\n"
    print("test_fractstep_discrete_laplacian example FAILED")


print(Text)
