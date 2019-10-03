{potflag}						##Record 1	
${nsurf0} ${nsurft} ${methflag} ${repflag}		##Record 2
${intflag}						##Record 3
if intflag == 0:		##Bulirsch-Stoer
    ${hstep0} ${eps} ${nprint}				##Record 3.0
elif intflag == 1:	##Rung-Kutta 4th order
    ${hstep} ${nprint}					##Record 3.1
${ranseed}						##Record 4
${ntraj} ${tflag1} ${tflag2} ${tflag3} ${tflag4}	##Record 5
if tflag1 == 1:
    ##nuclear momenta set to 0 every step
elif tflag1 == 2:					##Record 5.1
    ${ramptime} ${rampfact} ${nramp}	##temperature rescale
elif tflag1 == 3:					##Record 5.2
    ${andersen_temp} ${andersen_freq} ${scandth}  #Andersen thermostat

if tflag2 == 1:					##Record 5.3
    ${trajlist}

if tflag3 == 1:					##Record 5.4
    ${ntarget} ${ephoton} ${wphoton}
