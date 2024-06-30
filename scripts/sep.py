f = open("./CITCO_90_similarity_combined.sdf", "r")

currentfile = f.readline().strip() + ".sdf"
wf = open(currentfile, "a")
next = False

for line in f:
    if not next :
        if line.strip() != "$$$$":
            wf.write(line)
        else:
            next = True
    else:
        wf.close()
        currentfile = line.strip() + ".sdf"
        wf = open(currentfile, "a")
        wf.write(line)
        next = False
