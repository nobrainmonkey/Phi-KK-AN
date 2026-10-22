#!/bin/bash
# Build libUPCDst.so (dictionary + data classes of mc/upcdst) for local ROOT, no STAR env needed.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=$HERE/../mc/upcdst/StRoot/StUPCFilterMaker
OUT=$HERE/lib
mkdir -p "$OUT" && cd "$OUT"

CLASSES="StUPCEvent StUPCTrack StUPCBemcCluster StUPCVertex StUPCMCTrack StUPCMCEvent StUPCMCVertex
         StRPEvent StUPCRpsTrack StUPCRpsTrackPoint StUPCRpsCluster StUPCTofHit"
HDRS=""; SRCS=""
for c in $CLASSES; do HDRS="$HDRS $c.h"; SRCS="$SRCS $SRC/$c.cxx"; done

rootcling -f G__UPCDst.cxx -s libUPCDst.so -rml libUPCDst.so -rmf libUPCDst.rootmap \
  -I"$SRC" $HDRS "$HERE/../mc/upcdst/include/star-upcLinkDef.h"
c++ -shared -fPIC -O2 $(root-config --cflags) -I"$SRC" G__UPCDst.cxx $SRCS \
  $(root-config --libs) -lEG -o libUPCDst.so
echo "built $OUT/libUPCDst.so"
