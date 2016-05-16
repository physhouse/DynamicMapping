astyle --style=allman --indent=spaces=4 --pad-oper --pad-header --unpad-paren --align-pointer=name --align-reference=name --convert-tabs --recursive *.{cpp,h}

rm -f *.orig
rm -f */*.orig
rm -f */*/*.orig