${nmol} ${ezero}					##Record 6
for i in range(nmol):
    ${natom} ${initx} ${initp} ${initj} ${ezero_i}	##Record 7
    if initx == 0:					##Record 8.0
        for j in range(natom):
            ${target_mass_xyz} (${sym} ${mass} ${xx}
    elif initx == 1:					##Record 8.1
        ${rdum}						##Record 8.1.1
        for j in range(natom):
            ${target_mass} (${sym} ${mass})		##Record 8.1.2
    elif initx == 2:					##Record 8.2
        ${lreadhess} ${nmtype} ${nmqn}			##Record 8.2.1
        for j in range(natom):
            ${target_mass_xyz}				##Record 8.2.2
    elif initx == 3:					##Record 8.3
        for j in range(natom):
            ${target_mass}				##Record 8.3.1
        ${escatad} ${vvad} ${jjad} ${rrad} ${arrad}	##Record 8.3.2
    elif initx == 5:					##Record 8.5
        if temp0im >= 0:
            ${lreadhess} ${temp0im}			##Record 8.5.1
        elif temp0im < 0:
            ${lreadhess} ${temp0im} ${scale0im}		##Record 8.5.1
        for j in range(natom):
            ${target_mass_xyz}				##Record 8.5.3
    elif initx = 6:					##Record 8.6
        ${samptot} ${lbinsamp} ${sampfilexx} ${sampfilepp}	#Record 8.6.1
        for j in range(natom):
            ${target_mass}				##Record 8.6.3

    if initp == 0:					##Record 9.0
        ${temp0im} ${escale0im}

    if initj == 1:					##Record 10.1
        ${samptarg} ${sampjmin} ${sampjmax} ${sampjtemp1} ${samptjemp2} ${sampbrot1} ${sampbrot2}
        if samptarg < 0:				##Record 10.1.2
            ${ejsc}
