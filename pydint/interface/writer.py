"""
Executes the automation part of DiNT 
"""

import os
from mako.template import Template


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')

########## Run Info ##########

def dint_inp_1(potflag, nsurf0, nsurft, methflag, repflag, 
                  intflag, ranseed, ntraj,
                  tflag1, tflag2, tflag3, tflag4,
                  hstep0, eps, nprint,  # remaining are optional: enter as 0s if not used
                  hstep,
                  ramptime, rampfact, nramp,
                  andersen_temp, andersen_freq, scandth,
                  trajlist,
                  ntarget, ephoton, wphoton):
    """ writes the dint input file for trajectory info """
    # Set the dictionary for the DiNT input file
    fill_vals = {
        "potflag": potflag,
        "nsurf0": nsurf0,
        "nsurft": nsurft,
        "methflag": methflag,
        "repflag": repflag,
        "intflag": intflag,
        "tflag1": tflag1,
        "tflag2": tflag2,
        "tflag3": tflag3,
        "tflag4": tflag4
    }

    # Call other smaller functions if intflag == 0/1, tflag1 = 1/2, etc. and update dictionary
    fill_vals_intflag = dint_intflag_fill(fill_vals, intflag, hstep0, eps, nprint, hstep)
    fill_vals_tflag = dint_tflag_fill(fill_vals, tflag1, tflag2, tflag3, tflag4, 
                      ramptime, rampfact, nramp,
                      andersen_temp, andersen_freq, scandth,
                      trajlist, ntarget, ephoton, wphoton)

    fill_vals.update(fill_vals_intflag)
    fill_vals.update(fill_vals_tflag)

    # Set template name and path for the DiNT input file
    template_file_name = 'dint_inp_1.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build the DiNT input string
    input_str = Template(filename=template_file_path).render(**fill_vals)

    return input_str

def dint_intflag_fill(fillvals, intflag,    # remaining are optional: enter as 0s if not used
                        hstep0, eps, nprint,
                        hstep):
    """ writes the appropriate input flags for the given intflag value """

    if intflag==0:	# Bulirsch-Stoer integration
        fill_vals_intflag = {
            "hstep0": hstep0,
            "eps": eps,
            "nprint": nprint
        }
    elif intflag==1:	# Runge-Kutta 4th order integration
        fill_vals_intflag = {
            "hstep": hstep,
            "nprint": nprint
        }

    return fill_vals_intflag

def dint_tflag_fill(fillvals, tflag1, tflag2, tflag3, tflag4,
                      ramptime, rampfact, nramp,                # tflag1=2
                      andersen_temp, andersen_freq, scandth,    # tflag1=3
                      trajlist,                                 # tflag2=1
                      ntarget, ephoton, wphoton):               # tflag3=1
    """ writes the appropriate input flags for the given tflag values """

    if tflag1==2:	# temp rescale
        fill_vals_tflag = {
            "ramptime": ramptime,
            "rampfact": rampfact,
            "nramp": nramp
    }
    elif tflag1==3:	# Runge-Kutta 4th order integration
        fill_vals_tflag = {
            "andersen_temp": andersen_temp,
            "andersen_freq": andersen_freq,
            "scandth": scandth
    }

    if tflag2==1:
        fill_vals_tflag = {
            "trajlist": trajlist
        }

    if tflag3==1:
        fill_vals_tflag = {
            "ntarget": ntarget,
            "ephoton": ephoton,
            "wphoton": wphoton
        }
    
    return fill_vals_tflag

########## AG Info ##########

def dint_inp_2(nmol, ezero):
    """ writes the dint input file for trajectory info """
    # Set the dictionary for the DiNT input file
    fill_vals = {
        "nmol": nmol,
        "ezero": ezero
    }

    # Loop over vals for each AG and update dictionary
    for i in range(nmol):
        fill_vals_nmol = dint_inp_mol(fill_vals, natom, initx, initp, initj, ezero_i,
                sym, mass, xx, rdum,
                lreadhess, nmtype, nmqn,
                escatad, vvad, jjad, rrad, arrad,
                temp0im, scale0im,
                samptot, lbinsamp, sampfilexx, sampfilepp,
                samptarg, sampjmin, sampjmax,
                sampjtemp1, sampjtemp2, sampbrot1, sampbrot2, ejsc)

        fill_vals.update(fill_vals_nmol)

    # Set template name and path for the DiNT input file
    template_file_name = 'dint_inp_2.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build the DiNT input string
    input_str = Template(filename=template_file_path).render(**fill_vals)

    return input_str

def dint_inp_mol(fill_vals, natom, initx, initp, initj, ezero_i,
                sym, mass, xx, rdum,
                lreadhess, nmtype, nmqn,
                escatad, vvad, jjad, rrad, arrad,
                temp0im, scale0im,
                samptot, lbinsamp, sampfilexx, sampfilepp,
                samptarg, sampjmin, sampjmax,               # initj=1
                sampjtemp1, sampjtemp2, sampbrot1, sampbrot2, ejsc):
    """ loops over nmol times """
    fill_vals_nmol = {
        "natom": natom,
        "initx": initx,
        "initp": initp,
        "initj": initj,
        "ezero_i": ezero_i
    }

    if initx==0:
        fill_vals_x = {
            "sym": sym,
            "mass": mass,
            "xx": xx
        }
    elif initx==1:
        fill_vals_x = {
            "rdum": rdum,
            "sym": sym,
            "mass": mass
        }
    elif initx==2:
        fill_vals_x = {
            "lreadhess": lreadhess,
            "nmtype": nmtype,
            "nmqn": nmqn,
            "sym": sym,
            "mass": mass,
            "xx": xx
        }
    elif initx==3:
        fill_vals_x = {
            "sym": sym,
            "mass": mass,
            "escatad": escatad,
            "vvad": vvad,
            "jjad": jjad,
            "rrad": rrad,
            "arrad": arrad
        }
    elif initx==5:
        if temp0im >= 0:
            fill_vals_x = {
                "lreadhess": lreadhess,
                "temp0im": temp0im
            }
        elif temp0im < 0:
            fill_vals_x = {
                "lreadhess": lreadhess,
                "temp0im": temp0im,
                "scale0im": scale0im,
                "sym": sym,
                "mass": mass,
                "xx": xx
            }
    elif initx==6:
        fill_vals_x = {
            "samptot": samptot,
            "lbinsamp": lbinsamp,
            "sampfilexx": sampfilexx,
            "sampfilepp": sampfilepp,
            "sym": sym,
            "mass": mass
        }

    if initp==0:
        fill_vals_p = {
            "temp0im": temp0im,
            "scale0im": scale0im
        }

    if initj==1:
        fill_vals_j = {
            "samptarg": samptarg,
            "sampjmin": sampjmin,
            "sampjmax": sampjmax,
            "sampjtemp1": sampjtemp1,
            "sampjtemp2": sampjtemp2,
            "sampbrot1": sampbrot1,
            "sampbrot2": sampbrot2
        }
        if samptarg < 0:
            fill_vals_j.update({"ejsc": ejsc})

    fill_vals_nmol.update(fill_vals_x)
    fill_vals_nmol.update(fill_vals_p)
    fill_vals_nmol.update(fill_vals_j)

    return fill_vals_nmol

########## Collision Info ##########

def dint_inp_3(iorient, ldofrag, xx, pp, rel0qc, ttt, bminqc, bmaxqc, termflag, tnstep, tstime, tgradmag, tnoutcome, tsymb1, tsymb2, distcut, ioutput, ilist):
    """ writes the dint input file for trajectory info """
    # Set the dictionary for the DiNT input file
    fill_vals = {
        "iorient": iorient,
        "ldofrag": ldofrag,
        "termflag": termflag,
        "tnstep": tnstep,
        "ioutput": ioutput
    }

    # Call other smaller function if termflag == 1/2/3 and update dictionary
    fill_vals_iorient = dint_iorient_fill(fill_vals, iorient, xx, pp, rel0qc, ttt, bminqc, bmaxqc)  
    fill_vals_termflag = dint_termflag_fill(fill_vals, termflag, tstime, tgradmag, tnoutcome, tsymb1, tsymb2, distcut)
    
    fill_vals.update(fill_vals_iorient)
    fill_vals.update(fill_vals_termflag)

    #WIP
    if ioutput > 0:
        for i in range(ioutput):
            fill_vals_iout.update({"ilist"+"i": ilist[i]})

    # Set template name and path for the DiNT input file
    template_file_name = 'dint_inp_3.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build the DiNT input string
    input_str = Template(filename=template_file_path).render(**fill_vals)

    return input_str

def dint_iorient_fill(fillvals, iorient,
                        xx, pp,                         # iorient=0
                        rel0qc, ttt, bminqc, bmaxqc):   # iorient=1
    """ writes the appropriate input flags for the given tflag values """

    if iorient==0:
        fill_vals_iorient = {
            "xx": xx,
            "yy": yy
        }
    elif iorient==1:
        fill_vals_iorient = {
            "rel0qc": rel0qc,
            "ttt": ttt,
            "bminqc": bminqc,
            "bmaxqc": bmaxqc
        }
    
    return fill_vals_iorient

def dint_termflag_fill(fillvals, termflag,
                        tstime,                              # termflag=1
                        tgradmag,                            # termflag=2
                        tnoutcome, tsymb1, tsymb2, distcut): # termflag=3
    """ writes the appropriate input flags for the given tflag values """

    if termflag==1:
        fill_vals_termflag = {
            "tstime": tstime
        }
    elif termflag==2:
        fill_vals_termflag = {
            "tgradmag": tgradmag
        }
    elif termflag==3:
        fill_vals_termflag = {
            "tnoutcome": tnoutcome,
            "tsymb1": tsymb1,
            "tsymb2": tsymb2,
            "tsymb3": tsymb3,
            "distcut": distcut
        }
    
    return fill_vals_termflag

### WIP
def submission_script(drive_path, run_path, njobs):
    """ launches the job """

    # Write the bottom of the string
    job_exe_lines = '# Run several dint.x instances\n'
    job_exe_lines += 'cd {0}/run1\n'.format(run_path)
    job_exe_lines += 'time $DINTEXE < input.dat > output.dat &\n'
    for i in range(njobs-1):
        job_exe_lines += 'cd ../run{0}\n'.format(str(i+2))
        job_exe_lines += 'time $DINTEXE < input.dat > output.dat &\n'
    job_exe_lines += 'wait\n'

    # Set the dictionary for the DiNT input file
    fill_vals = {
        "job_exe_lines": job_exe_lines,
    }

    # Set template name and path for the DiNT input file
    template_file_name = 'submit.mako'
    template_file_path = os.path.join(drive_path, template_file_name)

    # Build the 1dmin input string
    sub_str = Template(filename=template_file_path).render(**fill_vals)

    return sub_str
