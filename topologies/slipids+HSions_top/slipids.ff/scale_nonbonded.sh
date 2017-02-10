#!/bin/bash
# script scales FF parameters for H-REM simulation in Gromacs

#This factor is actually gamma=gamma/gamma ; gamma=1/kT
GAMMA=0.31

awk -v gamma=$GAMMA '{
# if the line contains [*], put it into the variable brac
if ($0 ~ /\[.*\]/) brac = $0

# if the brac has word atomtypes inside and the line read doesnt contain the [*] header itself and doesnt contain comments (start with ;) and doesnt contain any other [*] header and does not start with # and if not blank (zero length)
if (brac ~ /\[ *atomtypes *\]/ && $0 !~ /\[ *atomtypes *\]/ && $0 !~ /^ *;/ && $0 !~ /^[ \t] *$/ && $0 !~ /^#/ && length($0) != 0 ) 
 #                                                                      charges           sigma-not, epsilon
 printf "%11s%5d%11.4f%9.3f%3s%13.5e%13.5e\n", $1"_sc",$2,$3,sqrt(gamma)*$4,$5,$6       , gamma*$7
else
 print $0

}' ffnonbonded.itp > ffnonbonded_scaled.itp

