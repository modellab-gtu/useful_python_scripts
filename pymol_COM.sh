#!/bin/bash
while read -r line
do

lig="$line"_out

#lig=DB8757_out

obabel "$lig".pdbqt -O "$lig".pdb -m

echo "100" > "$lig"_COM.xyz
echo "$lig" >> "$lig"_COM.xyz

for i in {1..100}
do

# Run PyMOL command and capture output
output=$(pymol -cq "$lig""$i".pdb -d 'select all; centerofmass sele')

# Parse output to extract x, y, and z coordinates
coordinates=$(echo "$output" | grep -oP '\[\s*(-?[0-9]+\.[0-9]+),\s*(-?[0-9]+\.[0-9]+),\s*(-?[0-9]+\.[0-9]+)\s*\]')
IFS=', ' read -r x y z <<<"${coordinates:1:${#coordinates}-2}"

# Print coordinates
#echo "Center of mass coordinates:"
echo "O" "$x" "$y" "$z" >> "$lig"_COM.xyz

done

done < list
