#add_library(ljlibs  #Lennard-Jones
#              LJ.f
#)

#add_library(tbpluslibs  #TB+exp/6
#              tbplusexp6.f
#              tbplusexp6dt.f
#              tbplusexp6all.f
#              tbplusexp6-2.f
#)

add_library(ch4libs #Methane
              ch4hALL.f
              ch4oALL.f
              ch4ohALL.f
              ch4o.f
              ch4oh.f
)

#install(TARGETS ch4libs DESTINATION ${CMAKE_LIBDIR})

#add_library(ch3flibs  #Fluoromethane
#              ch3f-rv2.f
#              ch3f-rv3aj2.f
#              ch3f-rv3aj2st.f
#              ch3f-rv3aj.f
#              ch3f-rv3.f
#              ch3f-rv4.f
#              ch3f-rv.f
#              ch3f-tmp.f
#)

#add_library(ho2libs  #HO2
#              ho2m.f
#              ho2h2oLS2.f
#              ho2heLS2.f
#              ho2heLS.f
#              ho2n2LS.f
#              ho2n2LS3C.f
#              ho2n2LSB.f
#              ho2n2LSC.f
#              ho2n2LSD.f
#              ho2n2LSE.f
#)

#add_library(fhlibs  #Li/NaFH
#              LiFHJ.f
#              LiFHJx.f
#              nafh2d.f
#              nafh2v.f
#              nafh2d-2d.f
#)

