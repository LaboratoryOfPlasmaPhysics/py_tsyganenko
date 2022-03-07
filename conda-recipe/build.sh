PATH=$PERFIX/bin:$PATH $BUILD_PREFIX/bin/python3 $(which meson) -Dbuildtype=release --prefix=$PREFIX  . build && cd build && ninja && ninja test && ninja install
