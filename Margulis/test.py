import subprocess

# values = (("r closed 1", "print geodesics 1"), ("r closed 2", "print geodesics 1"))

command = "snap"

mfld = 1
cutoff = 1.2

# lazy use of universal_newlines to prevent the need for encoding/decoding
p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, universal_newlines=True)  
output, err = p.communicate(input="r closed {}\n print geodesics {}\n".format(str(mfld), str(cutoff)))
# stderr is not connected to a pipe, so err is None
# we just want the result of the command
output = output.split('\n')
for item in output[2:-1]:
	print item[item.find("  ")+2:item.rfind("i")+1]

	# print(output[output.rfind(" "):-1])  # -1 to strip the final newline
