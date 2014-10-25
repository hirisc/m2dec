AC_INIT([m2dec],
        m4_esyscmd([./git-version-gen .tarball-version]),
        [bug@example])
LT_INIT
AM_INIT_AUTOMAKE(foreign dist-xz no-dist-gzip)
AC_CONFIG_SRCDIR(src/lib/mpeg2.cpp)
AC_ENABLE_SHARED
AC_DISABLE_STATIC
AC_PROG_LIBTOOL(libtool)
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AC_CONFIG_HEADER(src/lib/config.h)
dnl AC_PROG_LEX
dnl AC_PROG_YACC
dnl AC_CHECK_HEADERS(unistd.h)
AC_CHECK_SIZEOF(int *)
dnl AC_EXEEXT
AC_LANG_CPLUSPLUS
AC_C_BIGENDIAN
LDDISP=`sdl2-config --libs`
AC_SUBST([LDDISP])

AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debug info.],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=true])
AM_CONDITIONAL(NDEBUG, test x$debug = xfalse)

AC_ARG_ENABLE(fastdec,
[  --enable-fastdec    Turn on fast decoding mode.],
[case "${enableval}" in
  yes) fastdec=true ;;
  no)  fastdec=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-fastdec) ;;
esac],[fastdec=false])
AM_CONDITIONAL(FSTDEC, test x$fastdec = xtrue)

AC_ARG_ENABLE(display,
[  --enable-display    Turn on display using SDL.],
[case "${enableval}" in
  yes) display=true ;;
  no)  display=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-display) ;;
esac],[display=false])
AM_CONDITIONAL(DISP, test x$display = xtrue)

AC_ARG_ENABLE(display,
[  --enable-aalib    Turn on display using aalib.],
[case "${enableval}" in
  yes) aadisp=true ;;
  no)  aadisp=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-aalib) ;;
esac],[aadisp=false])
AM_CONDITIONAL(AALIB, test x$aadisp = xtrue)
dnl files to be generated
AC_OUTPUT(Makefile src/Makefile src/lib/Makefile src/app/Makefile Project/Makefile Project/VC/Makefile Project/HEW/Makefile Project/HEW/m2dec/Makefile Project/HEW/h264dec/Makefile Project/HEW/libm2dec/Makefile extras/Makefile)
