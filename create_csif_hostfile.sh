#!/bin/bash

# Output file
outfile="csif_hostfile"

# Original list of all CSIF computers
csif_pc=(pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 \
    pc11 pc12 pc13 pc14 pc15 pc16 pc17 pc18 pc19 pc20 \
    pc21 pc22 pc23 pc24 pc25 pc26 pc27 pc28 pc29 pc30 \
    pc31 pc32 pc33 pc34 pc35 pc36 pc37 pc38 pc39 pc40 \
    pc41 pc42 pc43)

# Delete output file
rm -f ${outfile}

# Try connecting to each computer
for pc in "${csif_pc[@]}";
do
    url="${pc}"
    echo -n "${pc}: "
    # Try to connect via SSH (should fail if computer is down or if NFS is down)
    nproc=$(ssh -o StrictHostKeyChecking=no \
        -o PreferredAuthentications=publickey "${url}" nproc)
    if [[ ! $? -eq 0 ]]; then
        # Failed
        echo "KO"
        continue
    fi

    # And to list of online computers
    echo "${url} slots=${nproc}" >> ${outfile}
    echo "OK"
done

# Shuffle list of online computers
cat ${outfile} | shuf -o ${outfile}

