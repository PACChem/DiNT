${iorient} ${ldofrag}				##Record 11
    if iorient == 0:				##Record 11.0
        for j in range(natom):
            ${xx}				##Record 11.0.1
        for j in range(natom):
            ${pp}				##Record 11.0.2
    elif iorient == 1:				##Record 11.1
        ${rel0qc} ${ttt} ${bminqc} ${bmaxqc}
${termflag} ${tnstep}				##Record 12
    if termflag == 1:				##Record 12.1
        ${tstime}
    elif termflag == 2:				##Record 12.2
        ${tgradmag}
    elif termflag == 3:				##Record 12.3
        ${tnoutcome}
        ${tsymb1} ${tsymb2} ${distcut}
    if ioutput > 0:
        ${ioutput} ${ilist1} ${ilist2} ..
