#/usr/bin/python
import Classes.Gadget2 as Gadget2
import sys

if len(sys.argv) != 3:
    print 'Proper use:\n python Gadget2ToAscii.py GadgetFileName OutputFile.txt'

G = Gadget2.ReadGadget2(sys.argv[1])#'1HqIso_Impact0_000' is the gadget2-output
G.WriteToAscii(sys.argv[2])



