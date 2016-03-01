#!/usr/bin/python
# coding: utf-8 -*-

import os
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt

def GCcount(sq,start,wind):
	at = 0
	sq += sq[:wind]
	end = start + wind
	for j in range(start,end):
		if(sq[j] == 'A' or sq[j] == 'T' or sq[j] == 'a' or sq[j] == 't'):
			at += 1.0
	return(at/wind)

def Acount(sq,start,wind):
	a = 0
	sq += sq[:wind]
	end = start + wind
	for j in range(start,end):
		if(sq[j] == 'A' or sq[j] == 'a'):
			a += 1.0
	return(a/wind)

def Arcount(sq,start,wind):
	a = 0
	at = 0
	sq += sq[:wind]
	end = start + wind
	for j in range(start,end):
		if(sq[j] == 'A' or sq[j] == 'a'):
			a += 1.0
		if(sq[j] == 'A' or sq[j] == 'T' or sq[j] == 'a' or sq[j] == 't'):
			at += 1.0
	return(a/at)
def main():
	fname = "Agrew_mito.fasta"
	i = nposl = 0
	window = 150
	step = 10
	'''
	window = ширина окна
	step = шаг. то есть длина графика будет во столько раз короче последовательности
	'''
	for seq_rec in SeqIO.parse(fname,"fasta"):
		nposl += 1
		dlina = seq_rec.seq.__len__()
		print dlina
		vals =list()
		avals = list()
		arvals = list()
		xvals = list()
		while(i<dlina):
			vals.append(GCcount(seq_rec.seq,i,window))
			avals.append(Acount(seq_rec.seq,i,window))
			arvals.append(Arcount(seq_rec.seq,i,window))
			xvals.append(i)
			i += step
		ngraph = len(vals) / 250 
		y0=min(vals)-.1
		y1=max(vals)+.1
		yy0 = list()
		yy1 = list()
		for j in range (0,250):
			yy0.append(.1)
			yy1.append(y1)
		plt.figure(1, figsize = (8, 12))
		for ii in range(0, ngraph):
			plt.subplot(ngraph,1,ii+1)
			x=xvals[ii*250:(ii+1)*250]
			y=vals[ii*250:(ii+1)*250]
			z=avals[ii*250:(ii+1)*250]
			t=arvals[ii*250:(ii+1)*250]
			plt.plot(x,y,'r-')
			plt.plot(x,z,'m-')
			plt.plot(x,t,'b.')
			plt.fill_between(x,yy1,y,alpha=.2,facecolor="Orange")
			plt.fill_between(x,yy0,z,alpha=.5,facecolor="Maroon")
			#plt.set_ylim=(y0,y1)
		plt.show()
		print y0,y1

	return()

if __name__=='__main__':
    main()
