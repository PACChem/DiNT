set(src_SOURCES 
 addrot.F 	andersen.F 	angmom.f 	ange.f 		atomdiatom.F 	bootstrap.f 
 bsstep.f 	carttojac.f 	checkfrus.f 	checkhop.F 	cn.f 		decocheck.F
 derivs.f 	detwell.f 	diamin.f 	diamin2.f 	diapot.f 	diapot2.f 
 diatom.F 	dint.F 		driver.f 	ejscale.f 	elecdeco.f 	ems.F
 ewkb.f 	ewkb2.f		finalstate.f 	getdvec2.f 	getgrad.f 	getpem.f 
 getrel.f 	getrho.f 	gettemp.f 	header.f 	hessems.F 	hop.f 
 honey.f 	initelec.f 	initmol.f 	lindemann.f 	mc.f 		mmid.f 
 noang.f 	normod.f 	nsphere.f 	nmturn.f 	nmpot.f 	orient.F 
 period.f 	period2.f 	pjsplit.f 	getplz.f 	popnorm.F 	preatomdiatom.f 
 prediatom.f 	premol.f 	preptraj.f 	pzextr.f 	radialdist.f 	ranclu.F 
 rancol.F 	rangaussian.F 	ranj.F 		ranno.F 	ranrot.F 	ransamp.F 
 ranstate.F 	rantherm.F 	ranwell.F 	rarray.f 	readin.f 	rk4.f 
 rtox.f 	spin.F 		stodeco.F 	takestep.f 	timing.f 	torphi.f 
 turn.f 	turn2.f 	vibwells.F 	vwkb.f 		vwkb2.f 	xptoy.f 
 airy.f
 ${PROJECT_SOURCE_DIR}/pot/${POT}.f
)

add_library(objs OBJECT ${src_SOURCES})
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/sprng/SRC")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/pot")
target_link_libraries(objs PUBLIC ${EXTRA_LIBS})
target_include_directories(objs PUBLIC "${PROJECT_BINARY_DIR}" ${EXTRA_INCLUDES})
install(TARGETS objs OBJECTS DESTINATION ${CMAKE_OBJDIR})
