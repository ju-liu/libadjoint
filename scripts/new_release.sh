# This is a stub for the steps required to create a new release

echo "Creating release tar ball"
PACKAGE=libadjoint
VERSION=2017.2.0
ARCHIVE=$PACKAGE-$VERSION.tar.gz
git archive --prefix=$PACKAGE-$VERSION/ -o $ARCHIVE $PACKAGE-$VERSION

