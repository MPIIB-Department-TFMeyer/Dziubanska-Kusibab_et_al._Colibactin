import fileinput

for line in fileinput.input():
    if line[0]=="#":
        print line.strip("\n\r")
        continue 
    fields = line.strip("\n\r").split("\t",-1)
    fields[0] = fields[0].replace("chr","")
    if fields[0]=="M":
        fields[0]="MT"
    print "\t".join(fields)


