${potflag}						##Record 1	
${nsurf0} ${nsurft} ${methflag} ${repflag}		##Record 2
${intflag}						##Record 3
if intflag == 0:		##Bulirsch-Stoer
    ${hstep0} ${eps} ${nprint}				##Record 3.0
elif intflag == 1:	##Rung-Kutta 4th order
    ${hstep} ${nprint}					##Record 3.1
${ranseed}						##Record 4
${ntraj} ${tflag1} ${tflag2} ${tflag3} ${tflag4}	##Record 5
if tflag1 == 0:
    ##nothing
elif tflag1 == 1:
    ##nuclear momenta set to 0 every step
elif tflag1 == 2:					##Record 5.1
    ${ramptime} ${rampfact} ${nramp}	##temperature rescale
elif tflag1 == 3:					##Record 5.2
    ${andersen_temp} ${andersen_freq} ${scandth}  #Andersen thermostat

if tflag2 == 0:
    ##nothing
elif tflag2 == 1:					##Record 5.3
    ${trajlist}

if tflag3 == 0:
    ##nothing
elif tflag3 == 1:					##Record 5.4
    ${ntarget} ${ephoton} ${wphoton}

if tflag4 == 0:
    ##nothing
elif tflag4 == 1:
    ##sd method

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

    if initp == -1:
	##see initx
    elif initp == 0:					##Record 9.0
        ${temp0im} ${escale0im}
    elif initp == 1:
        #initp = 0
    elif initp == 2:					##Record 9.2
        #read initp

    if initj == 0:
        #nothing
    elif initj == 1:					##Record 10.1
        ${samptarg} ${sampjmin} ${sampjmax} ${sampjtemp1} ${samptjemp2} ${sampbrot1} ${sampbrot2}
        if samptarg == 0:
            #do nothing
        elif samptarg < 0:				##Record 10.1.2
            ${ejsc}
        elif samptarg > 0:
            #scale vib momenta
##NMOL LOOP

${iorient} ${ldofrag}				##Record 11
    if iorient == 0:				##Record 11.0
        for j in range(natom):
            ${xx}				##Record 11.0.1
        for j in range(natom):
            ${pp}				##Record 11.0.2
    elif iorient == 1:				##Record 11.1
        ${rel0qc} ${ttt} ${bminqc} ${bmaxqc}
${termflag} ${tnstep}				##Record 12
    if termflag == 0:
        #run for ${tnstep}
    elif termflag == 1:				##Record 12.1
        ${tstime}
    elif termflag == 2:				##Record 12.2
        ${tgradmag}
    elif termflag == 3:				##Record 12.3
        ${tnoutcome}
        ${tsymb1} ${tsymb2} ${distcut}
    if ioutput == 0:				##Record 13
        ##write to all units
    elif ioutput > 0:
        ${ioutput} ${ilist1} ${ilist2} ..
