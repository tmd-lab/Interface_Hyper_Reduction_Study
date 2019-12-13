#!/bin/sh

if [ $# = 2 ]
then 
    echo "Correct call!"
    OUT=$2
elif [ $# = 1 ]
then
    echo "Acceptable call!"
    a="$1"
    OUT="${a%.*}.mat"
else
    echo "Wrong call - quitting!"
fi
echo "Preprocessing mtx files"
gawk 'BEGIN{start=0} ($1~/^\*M/){start=start+1;} ($1~/^\*\*/){start=0} (start==1 && $1!~/^\*/){print}' $1|gawk 'BEGIN{RS=",";ORS="\n"}{print}'|gawk '(NF!=0){print}' > ./.STIFFNESS.mtx
gawk 'BEGIN{start=0} ($1~/^\*M/){start=start+1;} ($1~/^\*\*/){start=0} (start==2 && $1!~/^\*/){print}' $1|gawk 'BEGIN{RS=",";ORS="\n"}{print}'|gawk '(NF!=0){print}' > ./.MASS.mtx
gawk 'BEGIN{start=0;op=0} ($1~/^\*M/){start=1} (start==1 && $1~/\*\*/){print $0}' $1|cut --complement -d ' ' -f1|gawk 'BEGIN{FS=".";OFS="."}(NF>=2){print}'|gawk 'BEGIN{RS=",\n";ORS=" "}{print}'|gawk '(NF>0){print}'|gawk 'BEGIN{RS=",";ORS=" "}{print}'|gawk '(NF>0){print}' > ./.RECOV.mtx
echo "Preprocessing mtx files done"

python <<EOF
import numpy as np
import scipy.io as io

Mv = np.loadtxt('.MASS.mtx');
Kv = np.loadtxt('.STIFFNESS.mtx');
R = np.loadtxt('.RECOV.mtx');

print("Reading mtx arrays done")

Nelm = len(Mv);
Nelk = len(Kv);
if (Nelm!=Nelk):
	sys.exit("GIGO - Mass & Stiffness not of same length.");
Nel = Nelm;

Nd = np.int((np.sqrt(1+8*Nel)-1)/2); # Solution of Nd(Nd+1)/2-Nel = 0

M = np.zeros((Nd,Nd));
K = np.zeros((Nd,Nd));

(xi,yi) = np.tril_indices(Nd);
M[xi,yi] = Mv;
M[yi,xi] = Mv;
K[xi,yi] = Kv;
K[yi,xi] = Kv;

print("Matrix extraction done - writing mat file")
dict = {"M": M, "K": K, "R": R};
io.savemat(".out.mat",dict);
print("Processing Over")
EOF
mv .out.mat $OUT
rm .STIFFNESS.mtx .MASS.mtx .RECOV.mtx
