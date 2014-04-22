
#!/bin/bash

# $1 = OnlyInteractions.txt
# example :10	1A6QA	0	2	145.372389
#10	1A6QA	5	6	136.443357
#10	1A6QA	7	8	122.130250
#10	1A6QA	7	9	-41.910383
#10	1A6QA	8	9	132.881214
#10	2E2RA	1	2	125.356397
#10	2E2RA	2	4	-77.816423
#10	2E2RA	0	5	92.877919
#10	2E2RA	2	5	144.185957
for i in 3 4 5 6 7 8 9 10
#for i in 3
do
	awk -v hlx=$i '{if($1 == hlx) print $2"\t"$3"\t"$4"\t"$5}' $1 > tmp
	python extractPatDiff.py tmp $i
	rm tmp
	echo $i
done
