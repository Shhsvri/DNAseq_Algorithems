#!/bin/python3

def z_array(s):
	assert len(s) > 1
	z = [len(s)] + [0] * (len(s)-1)
	for i in range(1, len(s)):
		
        
