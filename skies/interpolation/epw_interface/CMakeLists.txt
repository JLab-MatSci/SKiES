add_library(skies_epw_interface)

set_target_properties(skies_epw_interface PROPERTIES
		      EXPORT_NAME epw_interface
)

set_source_files_properties(epi_init.f90 epi_interpol.f90
			    PROPERTIES Fortran_PREPROCESS ON
)

target_sources(skies_epw_interface
PRIVATE
	epi_init.f90
	epi_interpol.f90
)

target_link_libraries(skies_epw_interface
PUBLIC
	qe_epw
	skies_quantities
#qe_mbd
	#qe_fftx
	#qe_fox
)

get_target_property(QE_PW_INCLUDE_DIR   qe_pw   INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(QE_MODULES_INCLUDE_DIR qe_modules INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(QE_UTILX_INCLUDE_DIR   qe_utilx   INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(QE_PHONON_PH_INCLUDE_DIR   qe_phonon_ph   INTERFACE_INCLUDE_DIRECTORIES)

target_include_directories(skies_epw_interface
PUBLIC
	${QE_PW_INCLUDE_DIR}
	${QE_MODULES_INCLUDE_DIR}
	${QE_UTILX_INCLUDE_DIR}
	${QE_PHONON_PH_INCLUDE_DIR}
)

target_link_directories(skies_epw_interface
PUBLIC
	 /usr/local/lib
)

#target_include_directories(skies_epw_interface
#PUBLIC
#	/home/ilyagaltsov/quantum_espresso/qe-7.1/build/Modules/mod/qe_modules
#)

#target_include_directories(skies_epw_interface
#PUBLIC
#		${QE_ROOT}/build/EPW/mod/qe_epw
#		${QE_ROOT}/build/Modules/mod/qe_modules
#		${QE_ROOT}/build/UtilXlib/mod/qe_utilx
#		${QE_ROOT}/build/PHonon/mod/qe_phonon_ph
#		${QE_ROOT}/build/PW/mod/qe_pw
#		${QE_ROOT}/build/external/fox/common/mod/FoX_fsys
#		${QE_ROOT}/build/external/fox/fsys/mod/FoX_fsys
#		${QE_ROOT}/build/external/fox/wxml/mod/FoX_fsys
#)

file(GLOB QE_LIBS "${QE_ROOT}/build/lib/*.a")
set(fftx_test_lib "${QE_ROOT}/build/lib/libqe_fftx_test_common.a")
list(FILTER QE_LIBS EXCLUDE REGEX "${fftx_test_lib}")

file(GLOB MDB_OBJS "${QE_ROOT}/build/external/mbd/src/CMakeFiles/mbd.dir/*.o")
file(GLOB FOX_DOM_OBJS "${QE_ROOT}/build/external/fox/dom/CMakeFiles/FoX_dom.dir/*.o")
file(GLOB FOX_FSYS_OBJS "${QE_ROOT}/build/external/fox/fsys/CMakeFiles/FoX_fsys.dir/*.o")
file(GLOB FOX_WXML_OBJS "${QE_ROOT}/build/external/fox/wxml/CMakeFiles/FoX_wxml.dir/*.o")
file(GLOB FOX_COMMON_OBJS "${QE_ROOT}/build/external/fox/common/CMakeFiles/FoX_common.dir/*.o")

#target_link_libraries(skies_epw_interface
#	PUBLIC
#		${QE_LIBS}
#		${MDB_OBJS}
#		${FOX_DOM_OBJS}
#		${FOX_FSYS_OBJS}
#		${FOX_WXML_OBJS}
#		${FOX_COMMON_OBJS}
#		${QE_ROOT}/build/FFTXlib/src/CMakeFiles/qe_fftx.dir/fft_ggen.f90.o
#	PUBLIC
#		${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

#if (SKIES_ENABLE_MPI)
#target_link_libraries(skies_epw_interface
#PUBLIC
#	MPI::MPI_C
#		#PRIVATE MPI::MPI_CXX
#	MPI::MPI_Fortran)
#endif()

install(TARGETS skies_epw_interface
	DESTINATION lib
	EXPORT skies-targets
)
