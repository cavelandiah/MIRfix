#!/usr/bin/env python

import unittest

from lib.Collection import *

class Test_Test_getfilename(unittest.TestCase):
    def test_true_0(self):
        self.assertEqual(getfilename("/usr/bin/perl.txt"), "perl")
    def test_true_1(self):
        self.assertEqual(getfilename("/usr/bin/perl"), "perl")
    def test_true_2(self):
        self.assertEqual(getfilename("perl"), "perl")
    def test_no_path(self):
        self.assertEqual(getfilename("421414"),"421414")
    def test_empty(self):
        self.assertLogs(getfilename(""))

class Test_Test_dofold(unittest.TestCase):
    # def dofold(listnewold,oldid,precseq,newid,newseq):
    def setUp(self):
        self.listnewold = []
        self.oldid = "Test1"
        self.precseq = "UGUGGGAUGAGGUAGUAGAUUGUAUAGUUUUAGGGUCAUACCCCAUCUUGGAGAUAACUAUACAGUCUACUGUCUUUCCCACG"
        self.newid = "Test2"
        self.newseq = "UGUGGGAUGAGGUAGUAGAUUGUAUAGUUUUAGGGUCAUACCCCAUCUUGGAGAUAACUAUACAGUCUACUGUCUUUCCCACG"
    
    def test_empty_sequences(self):
        with self.assertLogs(level="ERROR"):
            dofold(self.listnewold, self.oldid, "", self.newid, "")
    
    def test_valid_sequences(self):
        result = dofold(self.listnewold, self.oldid, self.precseq, self.newid, self.newseq)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 8)
        self.assertEqual(result[0], self.oldid)
        self.assertEqual(result[1], self.precseq)
        self.assertIsInstance(result[2], str)
        self.assertIsInstance(result[3], float)
        self.assertEqual(result[4], self.newid)
        self.assertEqual(result[5], self.newseq)
        self.assertIsInstance(result[6], str)
        self.assertIsInstance(result[7], float)

    def test_invalid_input(self):
        with self.assertLogs(level='ERROR'):
            dofold(self.listnewold, None, self.precseq, self.newid, self.newseq)

class Test_Test_foldnomat(unittest.TestCase):
    # def doalifold(alnfile,outdir):
    def test_foldnomat(self):
        inputfasta = "test_input.fasta"
        outputfasta = "test_output.fasta"

        with open(inputfasta, "w") as f:
            f.write(">test\n")
            f.write("UGUGGGAUGAGGUAGUAGAUUGUAUAGUUUUAGGGUCAUACCCCAUCUUGGAGAUAACUAUACAGUCUACUGUCUUUCCCACG")

        # Call the function
        foldnomat(inputfasta, outputfasta)
        # Read the output file
        with open(outputfasta, 'r') as f:
            output_contents = f.read()
            expected_output = ">test\nUGUGGGAUGAGGUAGUAGAUUGUAUAGUUUUAGGGUCAUACCCCAUCUUGGAGAUAACUAUACAGUCUACUGUCUUUCCCACG\n.((((((.((((((((((((((((((((((((((((........))))))))....)))))))))))))))))))))))))). (-39.00)"
            self.assertEqual(output_contents.strip(), expected_output)
        
        # Clean up test files
        os.remove(inputfasta)
        os.remove(outputfasta)

if __name__== 'main':
    unittest.main()
