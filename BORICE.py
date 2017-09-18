#!/usr/bin/env python2

"""
BORICE is software for the estimation of population mean outcrossing rates and
inbreeding coefficients using Bayesian methods.
Copyright (C) 2012 Vanessa A. Koelling, Patrick J. Monnahan, John K. Kelly

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
"""

import sys
import os
import math
import random
import csv
import time
from decimal import *

_seed_borice_rand = os.environ.get('BORICE_RAND_SEED')
if _seed_borice_rand:
	_seed_borice_rand = long(_seed_borice_rand) # the long function returns an integer of unlimited precision
	random.seed(_seed_borice_rand)
	
class Allele(object):
	"""An Allele is an object representing an allele in the population. It has a y value.
	"""
	def __init__(self, allele, locus):
		self.name = allele
		self.locus = locus
		self.y = 1.0
		self.af_list = []
	
	def __str__(self):
		"""Returns an allele string with the allele name and the locus name separated by a comma.
		"""
		return repr(self.name) + ', locus' + repr(self.locus)

class SingleLocusGenotype(object):
	"""A SingleLocusGenotype is an object made up of two alleles, 'first' and 'second'. Alleles in a genotype are ordered smallest (first) to largest (second).
	"""
	def __init__(self, first, second):
		self.first = min(first, second)
		self.second = max(first, second)
		self.imputed = False
		self.observed_imputed = False
		
	def __str__(self):
		"""Returns a genotype string with the first and second alleles separated by a slash.
		"""
		return repr(self.first) + '/' + repr(self.second)
	
	def calc_prob_offspring_given_selfing_mom_homozygote_standard_model(self, mom_g, locus):
		"""Calculates the probability of a homozygous offspring genotype given selfing and its maternal genotype; no null alleles, no allelic drop-out.
		"""
		if mom_g == None: # this is the case for a single-offspring family with no maternal genotype missing data at this locus; effectively skips the locus
			return 1.0
		else:
			mf = mom_g.first
			sf, ss = self.first, self.second
			sh = (sf == ss)
			if (sf == -9) and (ss == -9): # missing data; effectively skips the locus
				return 1.0
			else:
				if sh and (sf == mf): # offspring is homozygous
					return 1.0
				else: # impossible genotype
					return 0.0

	def calc_prob_offspring_given_selfing_mom_homozygote_null_model(self, mom_g, locus):
		"""Calculates the probability of a homozygous offspring genotype given selfing and its maternal genotype with null alleles.
		"""
		mf = mom_g.first
		sf, ss = self.first, self.second
		sh = (sf == ss)
		if (sf == -9) and (ss == -9): # missing data; effectively skips the locus
			return 1.0
		else:
			if (mf == 0): # mom is homozygous null
				return 0.0
			else: # mom is not homozygous null
				if sh: # offspring is homozygous
					if (sf == mf):
						return 1.0
					else: # impossible genotype
						return 0.0
				else: # offspring is het
					return 0.0
	
	def calc_prob_offspring_given_selfing_mom_heterozygote_standard_model(self, mom_g, locus):
		"""Calculates the probability of a heterozygous offspring genotype given selfing and its maternal genotype; no null alleles, no allelic drop-out.
		"""
		if mom_g == None: # this is the case for a single-offspring family with no maternal genotype missing data at this locus; effectively skips the locus
			return 1.0
		else:
			mf, ms = mom_g.first, mom_g.second
			sf, ss = self.first, self.second
			sh = (sf == ss)
			if (sf == -9) and (ss == -9): # missing data
				return 1.0
			elif (sf != mf) and (sf != ms) and (ss != mf) and (ss != ms): # impossible genotype
				return 0.0
			elif sh: # offspring is homozygous
				return 0.25
			elif (sf == mf) and (ss == ms): # offspring is identical het
				return 0.5
			else: # offspring is non-identical het
				return 0.0

	def calc_prob_offspring_given_selfing_mom_heterozygote_null_model(self, mom_g, locus):
		"""Calculates the probability of a heterozygous offspring genotype given selfing and its maternal genotype with null alleles.
		"""
		mf, ms = mom_g.first, mom_g.second
		sf, ss = self.first, self.second
		sh = (sf == ss)
		if (sf == -9) and (ss == -9): # missing data
			return 1.0
		else:
			if (mf == 0): # mom is het with null allele	
				if sh and (sf == ms): # offspring is homozygous and not null
					return 1.0
				else: # impossible genotype
					return 0.0
			else: # mom is het without null allele
				if (sf != mf) and (sf != ms) and (ss != mf) and (ss != ms): # impossible genotype
					return 0.0
				elif sh: # offspring is homozygous
					return 0.25
				elif (sf == mf) and (ss == ms): # offspring is identical het
					return 0.5
				else: # offspring is non-identical het
					return 0.0
	
	def calc_prob_offspring_given_outcrossing_mom_homozygote_standard_model(self, allele_list, allele_freq, mom_g, locus):
		"""Calculates the probability of a homozygous offspring genotype given outcrossing and its maternal genotype; no null alleles, no allelic drop-out.
		"""
		if mom_g == None: # this is the case for a single-offspring family with no maternal genotype missing data at this locus; effectively skips the locus
			return 1.0
		else:
			mf, ms = mom_g.first, mom_g.second
			sf, ss = self.first, self.second
		
			# missing data (-9) is not in the allele list and must be skipped
			if (sf == -9) and (ss == -9):
				return 1.0
			else:
				allele1 = allele_list.index(sf)
				allele2 = allele_list.index(ss)

			if (sf == mf):
				return allele_freq[allele2]
			elif (ss == mf):
				return allele_freq[allele1]
			else: # impossible genotype
				return 0.0

	def calc_prob_offspring_given_outcrossing_mom_homozygote_null_model(self, allele_list, allele_freq, mom_g, locus):
		"""Calculates the probability of a homozygous offspring genotype given outcrossing and its maternal genotype with null alleles.
		"""
		mf, ms = mom_g.first, mom_g.second
		sf, ss = self.first, self.second
		sh = (sf == ss)
		null_allele = allele_list.index(0)
		
		# missing data (-9) is not in the allele list and must be skipped
		if (sf == -9) and (ss == -9):
			return 1.0
		else:
			allele1 = allele_list.index(sf)
			allele2 = allele_list.index(ss)
		
		if (mf == 0): # mom is homozygous null
			if sh: # offspring is homozygous
				return (allele_freq[allele1])/(1.0 - allele_freq[null_allele])
			else: # offspring is het
				return 0.0
		else: # mom is homozygous not null
			if sh: # offspring is homozygous
				if (sf == mf):
					return (allele_freq[allele1] + allele_freq[null_allele])
				else: # impossible genotype
					return 0.0
			else: # offspring is het
				if (sf == mf):
					return allele_freq[allele2]
				elif (ss == mf):
					return allele_freq[allele1]
				else: # impossible genotype
					return 0.0
	
	def calc_prob_offspring_given_outcrossing_mom_heterozygote_standard_model(self, allele_list, allele_freq, mom_g, locus):
		"""Calculates the probability of a heterozygous offspring genotype given outcrossing and its maternal genotype; no null alleles, no allelic drop-out.
		"""
		if mom_g == None: # this is the case for a single-offspring family with no maternal genotype missing data at this locus; effectively skips the locus
			return 1.0
		else:
			mf, ms = mom_g.first, mom_g.second
			sf, ss = self.first, self.second
			sh = (sf == ss)
		
			# missing data (-9) is not in the allele list and must be skipped
			if (sf == -9) and (ss == -9):
				return 1.0
			else:
				allele1 = allele_list.index(sf)
				allele2 = allele_list.index(ss)
		
			if sh: # offspring is homozygote
				if (sf == mf):
					return allele_freq[allele1] * 0.5
				elif (sf == ms):
					return allele_freq[allele2] * 0.5
				else: # impossible genotype
					return 0.0
			else: # offspring is het
				if (sf != mf) and (sf != ms) and (ss != mf) and (ss != ms): # impossible genotype
					return 0.0
				elif (sf != mf) and (sf != ms):
					return allele_freq[allele1] * 0.5
				elif (ss != mf) and (ss != ms):
					return allele_freq[allele2] * 0.5
				else: # offspring is identical het
					return 0.5 * (allele_freq[allele1] + allele_freq[allele2])

	def calc_prob_offspring_given_outcrossing_mom_heterozygote_null_model(self, allele_list, allele_freq, mom_g, locus):
		"""Calculates the probability of a heterozygous offspring genotype given outcrossing and its maternal genotype with null alleles.
		"""
		mf, ms = mom_g.first, mom_g.second
		sf, ss = self.first, self.second
		sh = (sf == ss)
		null_allele = allele_list.index(0)
		
		# missing data (-9) is not in the allele list and must be skipped
		if (sf == -9) and (ss == -9):
			return 1.0
		else:
			allele1 = allele_list.index(sf)
			allele2 = allele_list.index(ss)
		
		if (mf == 0): # mom is het with null allele
			if sh:
				if (sf == ms): # offspring is homozygote and matches maternal allele 2
					return (allele_freq[allele1]/(1.0 - allele_freq[null_allele]) * 0.5) + ((allele_freq[allele1] + allele_freq[null_allele]) * 0.5)
				else: # offspring is homozygote and does not match maternal allele 2
					return allele_freq[allele1]/(1.0 - allele_freq[null_allele]) * 0.5
			else: # offspring is het
				if (sf == ms):
					return allele_freq[allele2] * 0.5
				elif (ss == ms):
					return allele_freq[allele1] * 0.5
				else: # impossible genotype
					return 0.0
		else: # mom is het without null allele
			if sh: # offspring is homozygous
				if (sf == mf):
					return (allele_freq[allele1] + allele_freq[null_allele]) * 0.5
				elif (sf == ms):
					return (allele_freq[allele1] + allele_freq[null_allele]) * 0.5
				else: # impossible genotype
					return 0.0
			else: # offspring is het
				if (sf != mf) and (sf != ms) and (ss != mf) and (ss != ms): # impossible genotype
					return 0.0
				elif (sf != mf) and (sf != ms):
					return allele_freq[allele1] * 0.5
				elif (ss != mf) and (ss != ms):
					return allele_freq[allele2] * 0.5
				else: # offspring is identical het
					return 0.5 * (allele_freq[allele1] + allele_freq[allele2])

	def calc_prob_mom(self, allele_list, allele_freq, inbreeding_coefficient):
		"""Calculates the probability of a maternal genotype given its inbreeding coefficient.
		"""
		sf, ss = self.first, self.second
		sh = (sf == ss)
		allele1 = allele_list.index(sf)
		allele2 = allele_list.index(ss)
		inb = (1.0 - inbreeding_coefficient)
		if sh:
			return (inb * math.pow(allele_freq[allele1], 2)) + (inbreeding_coefficient * allele_freq[allele1])
		else:
			return (inb * (2.0 * allele_freq[allele1] * allele_freq[allele2]))
				
	def impute_new_mom(self, allele_list, allele_freq, inbreeding_coefficient, locus, null_loci):
		"""Selects a new maternal genotype based on allele frequencies if there is no observed genotype at this locus. If the maternal genotype is an observed homozygote at this locus, a choice is made whether or not to step from the current genotype (homozygote or null heterozygote) based on the probability of those genotypes.
		"""
		assert round(sum(allele_freq), 1) == 1.0
		
		sf, ss = self.first, self.second
		sh = (sf == ss)
		
		if self.imputed:
			rand_num = random.random()
			cumulative_prob = 0
			i = 0
			null = null_loci[locus]
			if null:
				first_allele = 0
			else:
				first_allele = 1 # skip allele zero, which is at frequency 0.0
			# choose first allele
			while i == 0:
				cumulative_prob = cumulative_prob + allele_freq[first_allele]
				if rand_num < cumulative_prob:
					new_first = allele_list[first_allele]
					i = 1
				else:
					first_allele = first_allele + 1
		
			# choose second allele
			r_number = random.random()
			if r_number < inbreeding_coefficient:
				# mom is new homozygote
				new_second = new_first
			else:
				random_number = random.random()
				cum_prob = 0
				j = 0
				if null:
					second_allele = 0
				else:
					second_allele = 1 # skip allele zero, which is at frequency 0.0
				while j == 0:
					# mom is new het
					cum_prob = cum_prob + allele_freq[second_allele]
					if random_number < cum_prob:
						new_second = allele_list[second_allele]
						j = 1
					else:
						second_allele = second_allele + 1
		
			j = new_first
			k = new_second
			new_first = min(j, k)
			new_second = max(j, k)
		else:
			if self.observed_imputed: # this should only be true when the mom is an observed homozygote
				null_allele = allele_list.index(0)
				allele = allele_list.index(ss)
				r_num = random.random()
				inb = (1.0 - inbreeding_coefficient)
				if sh: # current genotype is same as observed homozygote
					null_het_prob = inb * (2.0 * allele_freq[allele] * allele_freq[null_allele]) # prob of null het
					if r_num < null_het_prob:
						new_first = 0
						new_second = ss
					else:
						new_first = sf
						new_second = ss
				else: # current genotype is null heterozygote alternative to observed homozygote
					homozygote_prob = (inb * math.pow(allele_freq[allele], 2)) + (inbreeding_coefficient * allele_freq[allele]) # prob of homozygote
					if r_num < homozygote_prob:
						new_first = ss
						new_second = ss
					else:
						new_first = sf
						new_second = ss
		return new_first, new_second

class Individual(object):
	"""An Individual object is a set of multilocus genotypes. Individuals belong to families, and are either an offspring or a mom of the family. Individuals also have inbreeding coefficients.
	"""
	def __init__(self, family, genotype_list, is_mom = False):
		self.family = family
		self.genotype_list = genotype_list
		self.inbreeding_coefficient = 0.0
		
		# designate individuals as either mom or offspring of a family
		if family is not None:
			if is_mom:
				family.add_mom(self)
			else:
				family.add_offspring(self)

	def __str__(self):
		"""Returns an individual's multilocus genotype as a string.
		"""
		x = "Genotype = ["
		for locus in self.genotype_list:
			x = x + str(locus) + ','
		return x + ']'
			
	def get_imputed_loci(self, population):
		"""Returns a list of genotypes that were imputed by def infer_mom or def tag_mom_genotype, and a list of the allele frequencies at each of those loci.
		"""
		imp = []
		alleles = []
		afreq = []
		for n, genotype in enumerate(self.genotype_list):
			allele_list = population.allele_list[n]
			allele_freq = population.allele_freq_list[n]
			if genotype == None:
				continue
			else:
				if (genotype.imputed == True) or (genotype.observed_imputed == True):
					imp.append(genotype)
					alleles.append(allele_list)
					afreq.append(allele_freq)
		assert len(imp) == len(alleles) == len(afreq)
		return imp, alleles, afreq
	
	def calc_inbreeding_coefficient(self, ih):
		"""Calculates the inbreeding coefficient (f) of an individual.
		"""
		assert (ih >= 0)
		if (ih >= 0) and (ih <= 5):
			self.inbreeding_coefficient = 1.0 - math.pow(0.5, ih)
		else:
			self.inbreeding_coefficient = 1.0
		return self.inbreeding_coefficient
	
	def calc_prob_offspring_geno(self, outcrossing_rate, population, mom, null_loci):
		"""Calculates an individual's multilocus genotype probability given its single-locus genotype probabilities.
		"""
		multilocus_selfing_prob = 1.0
		for n, genotype in enumerate(self.genotype_list):
			mom_g = mom.genotype_list[n]
			if mom_g == None:
				multilocus_selfing_prob = multilocus_selfing_prob * genotype.calc_prob_offspring_given_selfing_mom_homozygote_standard_model(mom_g, n)
			else:
				mf, ms = mom_g.first, mom_g.second
				mh = (mf == ms)
				null = null_loci[n]
				if null:
					if mh:
						multilocus_selfing_prob = multilocus_selfing_prob * genotype.calc_prob_offspring_given_selfing_mom_homozygote_null_model(mom_g, n)
					else:
						multilocus_selfing_prob = multilocus_selfing_prob * genotype.calc_prob_offspring_given_selfing_mom_heterozygote_null_model(mom_g, n)
				else:
					if mh:
						multilocus_selfing_prob = multilocus_selfing_prob * genotype.calc_prob_offspring_given_selfing_mom_homozygote_standard_model(mom_g, n)
					else:
						multilocus_selfing_prob = multilocus_selfing_prob * genotype.calc_prob_offspring_given_selfing_mom_heterozygote_standard_model(mom_g, n)
		
		#for testing only!
# 		for n, genotype in enumerate(self.genotype_list):
# 			mom_g = mom.genotype_list[n]
# 			mf, ms = mom_g.first, mom_g.second
# 			mh = (mf == ms)
# 			null = null_loci[n]
# 			if null:
# 				if mh:
# 					prob = genotype.calc_prob_offspring_given_selfing_mom_homozygote_null_model(mom_g, n)
# 				else:
# 					prob = genotype.calc_prob_offspring_given_selfing_mom_heterozygote_null_model(mom_g, n)
# 			else:
# 				if mh:
# 					prob = genotype.calc_prob_offspring_given_selfing_mom_homozygote_standard_model(mom_g, n)
# 				else:
# 					prob = genotype.calc_prob_offspring_given_selfing_mom_heterozygote_standard_model(mom_g, n)
# 			print("Genotype probability of selfing = %s" % prob)
# 			
		multilocus_outcrossing_prob = 1.0
		for n, genotype in enumerate(self.genotype_list):
			allele_list = population.allele_list[n]
			allele_freq = population.allele_freq_list[n]
			mom_g = mom.genotype_list[n]
			if mom_g == None:
				multilocus_outcrossing_prob = multilocus_outcrossing_prob * genotype.calc_prob_offspring_given_outcrossing_mom_homozygote_standard_model(allele_list, allele_freq, mom_g, n)
			else:
				mf, ms = mom_g.first, mom_g.second
				mh = (mf == ms)	
				null = null_loci[n]
				if null:
					if mh:
						multilocus_outcrossing_prob = multilocus_outcrossing_prob * genotype.calc_prob_offspring_given_outcrossing_mom_homozygote_null_model(allele_list, allele_freq, mom_g, n)
					else:
						multilocus_outcrossing_prob = multilocus_outcrossing_prob * genotype.calc_prob_offspring_given_outcrossing_mom_heterozygote_null_model(allele_list, allele_freq, mom_g, n)
				else:
					if mh:
						multilocus_outcrossing_prob = multilocus_outcrossing_prob * genotype.calc_prob_offspring_given_outcrossing_mom_homozygote_standard_model(allele_list, allele_freq, mom_g, n)
					else:
						multilocus_outcrossing_prob = multilocus_outcrossing_prob * genotype.calc_prob_offspring_given_outcrossing_mom_heterozygote_standard_model(allele_list, allele_freq, mom_g, n)
		
		#for testing only!
# 		for n, genotype in enumerate(self.genotype_list):
# 			allele_list = population.allele_list[n]
# 			allele_freq = population.allele_freq_list[n]
# 			mom_g = mom.genotype_list[n]
# 			mf, ms = mom_g.first, mom_g.second
# 			mh = (mf == ms)	
# 			null = null_loci[n]
# 			if null:
# 				if mh:
# 					prob = genotype.calc_prob_offspring_given_outcrossing_mom_homozygote_null_model(allele_list, allele_freq, mom_g, n)
# 				else:
# 					prob = genotype.calc_prob_offspring_given_outcrossing_mom_heterozygote_null_model(allele_list, allele_freq, mom_g, n)
# 			else:
# 				if mh:
# 					prob = genotype.calc_prob_offspring_given_outcrossing_mom_homozygote_standard_model(allele_list, allele_freq, mom_g, n)
# 				else:
# 					prob = genotype.calc_prob_offspring_given_outcrossing_mom_heterozygote_standard_model(allele_list, allele_freq, mom_g, n)
# 			print("Genotype probability of outcrossing = %s" % prob)
# 		
		selfing_rate = (1.0 - outcrossing_rate)
		prob_offspring_geno = (selfing_rate * multilocus_selfing_prob) + (outcrossing_rate * multilocus_outcrossing_prob)
		try:
			lnL = math.log(prob_offspring_geno)
		except:
			lnL = float('-inf')
		return lnL

	def calc_prob_mom_geno(self, population):
		"""Calculates a maternal individual's multilocus genotype probability given its single-locus genotype probabilities.
		"""
		multilocus_mom_prob = 1.0
		for n, genotype in enumerate(self.genotype_list):
			if genotype == None: # this is for the case where it is a single-offspring family with missing data and no maternal genotype
				multilocus_mom_prob = multilocus_mom_prob * 1.0
			else:
				allele_list = population.allele_list[n]
				allele_freq = population.allele_freq_list[n]
				multilocus_mom_prob = multilocus_mom_prob * genotype.calc_prob_mom(allele_list, allele_freq, self.inbreeding_coefficient)
		try:
			lnL = math.log(multilocus_mom_prob)
		except:
			lnL = float('-inf')
		return lnL

def tag_mom_genotype(momfirst, momsecond, offspring, locus_index, null_loci, family):
	"""Tags an observed maternal genotype as imputed if it is a homozygote, and returns a SingleLocusGenotype. This is for the purpose of dealing with null alleles.
	"""
# #	for testing only when moms need to be read in as is!
#  	slg = SingleLocusGenotype(momfirst, momsecond)
#  	slg.imputed = True
#  	return slg
	
	
	m_list = [momfirst, momsecond]
	works = True
	for child in offspring:
		cg = child.genotype_list[locus_index]
		# skips missing genotypes
		if (cg.first == -9) and (cg.second == -9):
			continue
		# checks that the observed genotype is possible based on the progeny genotype
		if (cg.first not in m_list) and (cg.second not in m_list):
			works = False
			break
	
	null = null_loci[locus_index]
	if null:
		if works:
			if (momfirst == momsecond):
				slg = SingleLocusGenotype(momfirst, momsecond)
				slg.observed_imputed = True
				return slg
			else:
				slg = SingleLocusGenotype(momfirst, momsecond)
				return slg
		else:
			if (momfirst == momsecond):
				null_allele = 0
				momfirst = null_allele
				m_list = [momfirst, momsecond]
				works = True
				for child in offspring:
					cg = child.genotype_list[locus_index]
					# skips missing genotypes
					if (cg.first == -9) and (cg.second == -9):
						continue
					# checks that imputed null genotype is possible based on progeny genotype
					if (cg.first == cg.second):
						if (cg.first not in m_list) and (null_allele not in m_list):
							works = False
							break
					else:
						if (cg.first not in m_list) and (cg.second not in m_list):
							works = False
							break
				if works:
					slg = SingleLocusGenotype(momfirst, momsecond)
					slg.observed_imputed = True
					return slg
				else:
					raise SingleLocusGenotypeError(cg.first, cg.second, locus_index, family)
			else:
				if (momfirst == 0):
					null_allele = 0
					m_list = [momfirst, momsecond]
					works = True
					for child in offspring:
						cg = child.genotype_list[locus_index]
						# skips missing genotypes
						if (cg.first == -9) and (cg.second == -9):
							continue
						# checks that null heterozygote genotype is possible based on progeny genotype
						if (cg.first == cg.second):
							if (cg.first not in m_list) and (null_allele not in m_list):
								works = False
								break
						else:
							if (cg.first not in m_list) and (cg.second not in m_list):
								works = False
								break
					if works:
						slg = SingleLocusGenotype(momfirst, momsecond)
						return slg
					else:
						raise SingleLocusGenotypeError(cg.first, cg.second, locus_index, family)
				else:
					raise SingleLocusGenotypeError(cg.first, cg.second, locus_index, family)
	else:
		if works:
			slg = SingleLocusGenotype(momfirst, momsecond)
			return slg
		else:
			raise SingleLocusGenotypeError(cg.first, cg.second, locus_index, family)

def find_mom_genotype(allele_set, offspring, locus_index, null_loci, family, valid_geno_index = 0):
	"""Imputes a maternal genotype, tags it as imputed, and returns a SingleLocusGenotype.
	"""
	# selects the first maternal genotype that works for the family
	null = null_loci[locus_index]
	if null:
		null_allele = 0
		for momfirst in allele_set:
			for momsecond in allele_set:
				m_list = [momfirst, momsecond]
				works = True
				for child in offspring:
					cg = child.genotype_list[locus_index]
					# skips missing genotypes
					if (cg.first == -9) and (cg.second == -9):
						continue
					# checks that imputed genotype is possible based on progeny genotype
					if (cg.first == cg.second):
						if (cg.first not in m_list) and (null_allele not in m_list):
							if valid_geno_index == 0:
								works = False
								break
							else:
								valid_geno_index = valid_geno_index - 1
					else:
						if (cg.first not in m_list) and (cg.second not in m_list):
							if valid_geno_index == 0:
								works = False
								break
							else:
								valid_geno_index = valid_geno_index - 1
				if works:
					slg = SingleLocusGenotype(momfirst, momsecond)
					slg.imputed = True
					return slg
		
		if works:
			pass
		else:
			raise SingleLocusGenotypeError(cg.first, cg.second, locus_index, family)
	else:
		for momfirst in allele_set:
			for momsecond in allele_set:
				m_list = [momfirst, momsecond]
				works = True
				for child in offspring:
					cg = child.genotype_list[locus_index]
					# skips missing genotypes
					if (cg.first == -9) and (cg.second == -9):
						continue
					# checks that imputed genotype is possible based on progeny genotype
					if (cg.first not in m_list) and (cg.second not in m_list):
						if valid_geno_index == 0:
							works = False
							break
						else:
							valid_geno_index = valid_geno_index - 1
				if works:
					slg = SingleLocusGenotype(momfirst, momsecond)
					slg.imputed = True
					return slg
		if works:
			pass
		else:
			raise SingleLocusGenotypeError(cg.first, cg.second, locus_index, family)

class Family(object):
	"""A Family object is a maternal individual, its offspring, and its inbreeding history.
	"""
	def __init__(self, name):
		self.name = name
		self.mom = None
		self.pop_name = None
		self.offspring = []
		self.inbreeding_history_list = []
		self.inbreeding_history = 0
		self.locus_genotypes = []
		self.possible_genotypes = []
	
	def add_mom(self, mom):
		"""Adds a mom to a family unless one's already there.
		"""
		assert self.mom is None
		self.mom = mom
				
	def add_offspring(self, offspring):
		"""Adds offspring to a family."""
		self.offspring.append(offspring)
	
	def infer_mom(self, null_loci):
		"""Infers a maternal genotype for a family from offspring data.
		"""
		# constructs a list of observed alleles
		first_child = self.offspring[0]
		num_loci = len(first_child.genotype_list)
		observed_alleles = []
		for i in range(num_loci):
			observed_alleles.append(set())

		# adds observed alleles at each locus 
		for child in self.offspring:
			for i in range(num_loci):
				cg = child.genotype_list[i]
				allele_set = observed_alleles[i]
				allele_set.add(cg.first)
				allele_set.add(cg.second)
				allele_set.discard(-9)
				null = null_loci[i]
				if null:
					allele_set.add(0)

		# case for no mom genotype
		if not self.mom:
			mom_geno_list = []
			for i in range(num_loci):
				allele_set = observed_alleles[i]
				if len(allele_set) == 0:
					mg = None
				else:
					mg = find_mom_genotype(allele_set, self.offspring, i, null_loci, self.name)
					assert mg
				mom_geno_list.append(mg)
			assert len(mom_geno_list) == num_loci
			the_family = self
			self.mom = Individual(the_family, mom_geno_list, is_mom = True)
		# case for partial mom genotype; some loci inferred
		else:
			mgl = self.mom.genotype_list
			missing = []
			for n, geno in enumerate(mgl):
				if (geno.first == -9) and (geno.second == -9):
					missing.append(n)
				else: # case for observed mom, but null allele possible
					mg = tag_mom_genotype(geno.first, geno.second, self.offspring, n, null_loci, self.name)
					assert mg
					self.mom.genotype_list[n] = mg
					
			for i in missing: # this section imputes a possible genotype for the missing loci
				allele_set = observed_alleles[i]
				if len(allele_set) == 0:
					mg = None
				else:
					mg = find_mom_genotype(allele_set, self.offspring, i, null_loci, self.name)
					assert mg
				self.mom.genotype_list[i] = mg
				
	def __str__(self):
		"""Returns a string identifying a family and its maternal individual.
		"""
		r = ["\nFamily = %s \nMaternal %s" % (self.name, self.mom)]
		for o in self.offspring:
			r.append(str(o))
		return "\nProgeny ".join(r)

	def calc_mom_lnL(self):
		"""Calculates the ln likelihood value for the mom of a family.
		"""
		lnL = self.mom.calc_prob_mom_geno(self.population_name)
		return lnL
		
	def calc_family_lnL(self, outcrossing_rate, null_loci):
		"""Calculates the ln likelihood value for the entire family (mom and offspring).
		"""
		lnL = 0.0
		for offspring in self.offspring:
			lnL = lnL + offspring.calc_prob_offspring_geno(outcrossing_rate, self.population_name, self.mom, null_loci)
		lnL = lnL + self.mom.calc_prob_mom_geno(self.population_name)
		return lnL
	
	def calc_progeny_lnL(self, outcrossing_rate, null_loci):
		"""Calculates the ln likelihood value for only the offspring of a family.
		"""
		lnL = 0.0
		for offspring in self.offspring:
			lnL = lnL + offspring.calc_prob_offspring_geno(outcrossing_rate, self.population_name, self.mom, null_loci)
		return lnL
		
class Population(object):
	"""A Population object consists of multiple families. Each population also has a list of allele frequencies specific to that population, as well as an outcrossing rate.
	"""
	def __init__(self, allele_list, allele_freq_list, y_values, outcrossing_rate):
		self.allele_list = allele_list
		self.allele_freq_list = allele_freq_list
		self.family_list = []
		self.outcrossing_rate = float(outcrossing_rate)
		self.ih_prob_list = []
		self.y_values = y_values
	
	def calculate_new_af_list(self, locus_alleles, locus_index, step, burn_in, null_loci):
		"""Calculates a new list of allele frequencies at a locus based on the y-value of each allele.
		"""
		null = null_loci[locus_index]
		if null:
 			for n, allele in enumerate(locus_alleles):
 				self.y_values[locus_index].append(allele.y)
 			
 			new_af_list = []
 			for n, allele in enumerate(locus_alleles):
 				new_af = allele.y / sum(self.y_values[locus_index])
 				new_af_list.append(new_af)
 				if step > burn_in:
 					if (step % 10) == 0:
 						allele.af_list.append(new_af)
 		else:
			for n, allele in enumerate(locus_alleles):
				if n == 0:
					continue
				else:
					self.y_values[locus_index].append(allele.y)
			
			new_af_list = [0.0]
			for n, allele in enumerate(locus_alleles):
				if n == 0:
					continue
				else:
					new_af = allele.y / sum(self.y_values[locus_index])
					new_af_list.append(new_af)
					if step > burn_in:
						if (step % 10) == 0:
							allele.af_list.append(new_af)
		return new_af_list
	
	def calc_pop_lnL(self, null_loci):
		"""Calculates the ln likelihood of a population summed over families.
		"""
		lnL = 0.0
		for family in self.family_list:
			lnL = lnL + family.calc_family_lnL(self.outcrossing_rate, null_loci)
		return lnL
		
	def calc_ih_prob(self):
		"""Calculates a list of inbreeding history probabilities based on the population outcrossing rate.
		"""
		outcrossing_rate = self.outcrossing_rate
		for ih in range(0,7):
			if ih == 0:
				ih_prob = outcrossing_rate
				self.ih_prob_list.append(ih_prob)
			else:
				if (ih > 0) and (ih < 6):
					ih_prob = (math.pow(1 - outcrossing_rate, ih)) * outcrossing_rate
					self.ih_prob_list.append(ih_prob)
				else:
					ih_prob = 1 - sum(self.ih_prob_list)
					self.ih_prob_list.append(ih_prob)

class SingleLocusGenotypeError(Exception):
	"""Makes a SingleLocusGenotypeError class. If an impossible genotype is encountered in the data, it prints the genotype and identifies the family in which it occurs.
	"""
	def __init__(self, first, second, locus, family):
		self.first = first
		self.second = second
		self.locus_number = locus
		self.family = family
		
	def __str__(self):
		return "Impossible genotype %s/%s at locus %s in family %s!" % (self.first, self.second, self.locus_number + 1, self.family)

class CSVFileParseException(Exception):
	"""Makes a CSVFileParseException class.
	"""
	def __init__(self, stream, line, msg):
		self.filename = stream.name
		self.line = line
		self.msg = msg
	
	def __str__(self):
		return "Error parsing %s at line %d:\n	%s\n" % (self.filename, self.line, self.msg)

# parsing function for csv file containing genotype data for families and populations
def parse_csv(stream, sep):
	"""Reads a CSV file and returns marker names, families, and genotypes.
	"""
	import csv
	reader = csv.reader(stream, delimiter = sep) # file object passed in
	line_iterator = iter(reader)
	
	# first line should be number of markers, indicator of presence of population name, indicator of presence of a subgroup name
	first_row = line_iterator.next()
	if len(first_row) < 3:
		raise CSVFileParseException(stream, 1, "Expecting at least three columns in the first row")

	marker_number = first_row[0]
	try:
		num_markers = int(marker_number)
	except:
		raise CSVFileParseException(stream, 1, 'Expecting the number of markers in the first cell of the first row (Found "%s")' % marker_number)

	pop_indicator = first_row[1]
	if pop_indicator == '1':
		pop_name = True
	elif pop_indicator == '0':
		pop_name = False
	else:
		raise CSVFileParseException(stream, 1, 'Expecting a 0 or 1 in the second column of the first row (Found "%s") to indicate population names absent or present' % pop_indicator)

	subgroup_indicator = first_row[2]
	if subgroup_indicator == '1':
		subgroup_name = True
		raise CSVFileParseException(stream, 1, 'Subgroup names not supported')
	elif subgroup_indicator == '0':
		subgroup_name = False
	else:
		raise CSVFileParseException(stream, 1, 'Expecting a 0 or 1 in the third column of the first row (Found "%s") to indicate subgroup names absent or present' % subgroup_indicator)

	second_row = line_iterator.next()
	if len(second_row) < num_markers:
		raise CSVFileParseException(stream, 2, "Expecting at least %s columns of marker names in the second row")
	# the lists within marker_names each contain two index positions, index 0 = marker name, index 1 = empty
	marker_names = []
	for n, marker in enumerate(second_row):
		if n >= num_markers:
			break
		if not marker.strip():
			raise CSVFileParseException(stream, 2, "Found an empty cell in column %d of line 2 (expected a marker name)" % (1 + n))
		marker_names.append([marker.strip(), []]) # add the marker names from the second row of the file to the list of marker names
	
	
	expected_column_number = 2 + 2*num_markers
	families_in_pop = {}
	for n, row in enumerate(line_iterator):
		if not row:
			continue # allow blank lines by skipping them
		if len(row) < expected_column_number:
			raise CSVFileParseException(stream, 3 + n, "Expecting at least %d columns, but found only" % (expected_column_number, len(row)))
		
		# the below code tags rows as containing either mom or offspring individuals
		population_name = row[1] # the population name is in the second cell of the row
		family_name = row[0] # the family name is in the first cell of the row
		if family_name.endswith('!'): # moms are indicated with an ! at the end of the family name
			mom = True
			family_name = family_name[0:-1] # the mom's family name is the text minus the !
		else:
			mom = False

		# puts family names in the name_to_fam dictionary
		# binds families to their population; family names (such as numeric values) can be reused without error
		key = (population_name, family_name)
		if key in families_in_pop:
			family = families_in_pop[key]
		else:
			family = Family(family_name)
			family.pop_name = population_name
			families_in_pop[key] = family

		assert population_name == family.pop_name

		offset = 2	# the first 2 columns in the csv file are not genotype data
		genotype_list = []
		for i in range(num_markers):
			first = None
			second = None
			
			# finds the first allele of a locus
			try:
				first = int(row[offset])
			except:
				if row[offset] != '?':
					raise CSVFileParseException(stream, 3 + n, "Expecting a number for an allele in column %d, but found %s" % (offset + 1, row[offset]))
			
			# finds the second allele of a locus
			try:
				second = int(row[offset + 1])
			except:
				if row[offset + 1] != '?':
					raise CSVFileParseException(stream, 3 + n, "Expecting a number for an allele in column %d, but found %s" % (offset + 2, row[offset + 1]))
			
			# constructs a multilocus genotype_list
			if (first is not None) and (second is not None):
				genotype_list.append(SingleLocusGenotype(first, second))
			else:
				genotype_list.append(None)
			alleles_per_locus = [first, second]
			info_for_this_locus = marker_names[i]
			info_for_this_locus[1] = (alleles_per_locus + info_for_this_locus[1])
			offset = offset + 2
		individual = Individual(family, genotype_list, mom)
	return marker_names, families_in_pop.values()	

class BORICE(object):
	"""An object that encompasses the main functional components of the BORICE software.
	"""
	def __init__(self):
		#Progress of calculation
		self.current_step = 0		

	def getStep(self):
		return self.current_step

	def main(self, file_name, locus_model, num_steps = None, burn_in = None, outcrossing_rate_tuning_parameter = 0.05, allele_freq_tuning_parameter = 0.1, initial_outcrossing_rate = 0.5, writeOutput2 = True, writeOutput3 = True, writeOutput4 = True):
		input = open(file_name, 'rU')
		#print writeOutput2, writeOutput3, writeOutput4
		try:
			marker_names, families = parse_csv(input, ',')	
		except CSVFileParseException as x:
			sys.exit(str(x))
		
		assert len(marker_names) == len(locus_model)
		# main body of code begins here
		sorted_alleles_all_loci = []
		allele_freq_all_loci = []
		initial_y_values_all_loci = []
		index = 0
		for marker in marker_names:
			null = locus_model[index]
			# alleles at each locus are identified and missing values eliminated
			num_alleles = marker[1]
			unique_alleles = set(num_alleles)
			unique_alleles.discard(-9)
			#unique_alleles.discard(0)
			# allele list is created and sorted
			sorted_alleles = [0] # need to start the list with zero here; 0 is the null allele (not observed)
			al = list(sorted(unique_alleles))
			sorted_alleles.extend(al)
			sorted_alleles_all_loci.append(sorted_alleles)
			# initial allele frequencies are set
			if null:
				alleles = len(sorted_alleles)
 				freq = 1.0 / alleles
 				assert freq
 				af = [] # the initial allele frequency list; allele zero is the null allele
 			else:
				alleles = len(unique_alleles)
				freq = 1.0 / alleles
				assert freq
				af = [0.0] # if no null allele, then allele zero has a frequency of zero
			for i in range(alleles):
				af.append(freq)
			assert af
			allele_freq_all_loci.append(af)
			# makes a list of initial y values; initial y is 1.0 for each allele
			y_values = []
			initial_y_values_all_loci.append(y_values)
			index = index + 1
		assert len(sorted_alleles_all_loci) == len(allele_freq_all_loci) == len(initial_y_values_all_loci)
		population = Population(sorted_alleles_all_loci, allele_freq_all_loci, initial_y_values_all_loci, initial_outcrossing_rate)
		
		#print(population.allele_list)
		#print(population.allele_freq_list)
		#print(population.y_values)
		
		all_alleles = []
		for n, locus in enumerate(population.allele_list):
			alleles = population.allele_list[n]
			locus_alleles = []
			for i, allele in enumerate(alleles):
				allele_object = Allele(allele, locus)
				locus_alleles.append(allele_object)
			all_alleles.append(locus_alleles)

		families = sorted(families)
		# calls the function to initially impute maternal genotypes; also lists families in a population
		for fam in families:
			fam.infer_mom(locus_model) #this gives the initial inference of missing maternal loci; tags loci as imputed
			fam.population_name = population
			population.family_list.append(fam)
			#print(fam)
			#mom_lnL = fam.calc_mom_lnL()
			#print("Log-likelihood of mom = %s" % mom_lnL)
			#progeny_lnL = fam.calc_progeny_lnL(population.outcrossing_rate, locus_model)
			#print("Log-likelihood of progeny = %s" % progeny_lnL)
			#print("Family = %s mom likelihood = %s progeny likelihood = %s" % (fam.name, mom_lnL, progeny_lnL))
			# the below function is called here to catch errors in the progeny genotypes prior to imputing parameters
			#for child in fam.offspring:
			#	child.calc_prob_offspring_geno(fam.population_name.outcrossing_rate, fam.population_name, fam.mom, is_stepping)
		
		#print(population.allele_freq_list)
		#prev_lnL = population.calc_pop_lnL(locus_model)
		#print("Log-likelihood of data = %s" % prev_lnL)
		
		# creates a list of lists for each family for storage of imputed genotypes at each locus
		for fam in families:
			for n, locus in enumerate(fam.mom.genotype_list):
				fam.locus_genotypes.append([])
				fam.possible_genotypes.append([])

		#creates four output files
		borice_output1 = open('BORICE_output1.txt', 'wb')
		if(writeOutput2):
			borice_output2 = open('BORICE_output2.txt', 'wb')
		if(writeOutput3):
			borice_output3 = open('BORICE_output3.txt', 'wb')
			borice_output3.write("List of t, F, and ln likelihoood values from every 10 steps in the chain beyond the burn-in\nt\tF\tLn Likelihood of the Data\n")
		if(writeOutput4):
			borice_output4 = open('BORICE_output4.txt', 'wb')
	
		# below lists needed for storage of t, ih, and F before output to text files
		t_list = []
		ih_list = []
		F_list = []
		pop_lnL_list = []
		
		start_time = time.asctime()
		print("start time was %s" % start_time)
		
		# below is the code to perform Bayesian inference of outcrossing rate (t), inbreeding history, allele frequencies, maternal genotypes
		if burn_in is None:		
			burn_in = 999
		if num_steps is None:		
			num_steps = 100000
                #print num_steps
		for step in range(num_steps):
			if step != 0:
				self.current_step = step - 1
                        #print self.current_step
			
			prev_t = population.outcrossing_rate
			prev_lnL = population.calc_pop_lnL(locus_model)
			#print(prev_lnL)
			# changes outcrossing rate
			t_prime = (prev_t + ((random.random() - 0.5) * float(outcrossing_rate_tuning_parameter)))
			if t_prime < 0.0:
				t_prime = (0.0 - t_prime)
			if t_prime > 1.0:
				t_prime = (2.0 - t_prime)
			population.outcrossing_rate = t_prime
			#print(population.outcrossing_rate)
			lnL = population.calc_pop_lnL(locus_model)
			#print(lnL)
			if (lnL == float('-inf')):
				population.outcrossing_rate = prev_t
				prev_lnL = prev_lnL
				#print("1")
			else:
				lnL_ratio = (lnL - prev_lnL)
				if (lnL_ratio > 0):
					prev_t = population.outcrossing_rate
					prev_lnL = lnL
					#print("2")
				else:
					random_number = random.random()
					value = math.exp(lnL_ratio)
					if (random_number < value):
						prev_t = population.outcrossing_rate
						prev_lnL = lnL
						#print("2")
					else:
						population.outcrossing_rate = prev_t
						prev_lnL = prev_lnL
						#print("1")
	
			if step > burn_in:
				if (step % 10) == 0:
					t_list.append(population.outcrossing_rate)
					if(writeOutput3):
						borice_output3.write("%.2f" % population.outcrossing_rate + "\t")
			
			# changes inbreeding history
			population.calc_ih_prob()
			f_list = []
			for fam in families:
				#print(fam)
				prev_ih = fam.inbreeding_history
				#print(prev_ih)
				prev_f = fam.mom.calc_inbreeding_coefficient(fam.inbreeding_history)
				#print(prev_f)
				prev_mom_lnL = fam.mom.calc_prob_mom_geno(fam.population_name)
 				#print(prev_mom_lnL)
 				
				rand_num = random.random()
				cumulative_prob = 0
				ih_value = 0
				i = 0
				# choose new inbreeding history value
				while i == 0:
					cumulative_prob = cumulative_prob + population.ih_prob_list[ih_value]
					if rand_num < cumulative_prob:
						new_ih = ih_value
						i = 1
					else:
						ih_value = ih_value + 1
 		
				fam.inbreeding_history = new_ih
				#print(fam.inbreeding_history)
				new_f = fam.mom.calc_inbreeding_coefficient(fam.inbreeding_history)
				#print(new_f)
				lnL = fam.mom.calc_prob_mom_geno(fam.population_name)
				#print(lnL)
 				if (lnL == float('-inf')):
 					fam.inbreeding_history = prev_ih
					prev_f = fam.mom.calc_inbreeding_coefficient(fam.inbreeding_history)
					f_list.append(prev_f)
					#print("1")
				else:
					lnL_ratio = (lnL - prev_mom_lnL)
					if (lnL_ratio > 0):
						prev_ih = fam.inbreeding_history
						new_f = fam.mom.calc_inbreeding_coefficient(fam.inbreeding_history)
						f_list.append(new_f)
						#print("2")
					else:
						rand = random.random()
						value = math.exp(lnL_ratio)
						if (rand < value):
							prev_ih = fam.inbreeding_history
							new_f = fam.mom.calc_inbreeding_coefficient(fam.inbreeding_history)
							f_list.append(new_f)
							#print("2")
						else:
							fam.inbreeding_history = prev_ih
							prev_f = fam.mom.calc_inbreeding_coefficient(fam.inbreeding_history)
							f_list.append(prev_f)
							#print("1")
 		
				if step > burn_in:
					if (step % 10) == 0:
						ih_list.append(fam.inbreeding_history)
						fam.inbreeding_history_list.append(fam.inbreeding_history)
 			
			pop_inbreeding_coefficient = sum(f_list)/len(f_list)

			if step > burn_in:
				if (step % 10) == 0:
					F_list.append(pop_inbreeding_coefficient)
					if(writeOutput3):
						borice_output3.write("%.2f" % pop_inbreeding_coefficient + "\t")

			del population.ih_prob_list[:]
 			
 			if (step % 10) == 0:
				# changes allele frequencies
				for locus_index, locus in enumerate(all_alleles):
					allele_freq = population.allele_freq_list[locus_index]
					prev_allele_freq = list(allele_freq)
					#print(prev_allele_freq)
					locus_alleles = all_alleles[locus_index]
					null = locus_model[locus_index]
					if null:
						random_allele = random.randint(0, len(locus_alleles) - 1)
						allele = locus_alleles[random_allele]
					elif len(allele_freq) == 1: # this is to account for single-offspring families with missing data at a locus with no null alleles
						continue
					else:
						random_allele = random.randint(1, len(locus_alleles) - 1) # skips allele zero, which should remain at a frequency of zero with no null alleles
						allele = locus_alleles[random_allele]
						#print(allele)
						
					prev_y = allele.y
					#print(prev_y)
					prev_lnL = population.calc_pop_lnL(locus_model)
					#print(prev_lnL)
					new_y = (prev_y + ((random.random() - 0.5) * float(allele_freq_tuning_parameter)))
					if new_y < 0:
						new_y = (0.0 - new_y)
					allele.y = new_y
					#print(allele.y)
					new_af_list = population.calculate_new_af_list(locus_alleles, locus_index, step, burn_in, locus_model)
					#print(new_af_list)
					#calculates new lnL based on new allele frequencies
					population.allele_freq_list[locus_index] = new_af_list	
					lnL = population.calc_pop_lnL(locus_model)
					#print(lnL)
					population.y_values[locus_index] = []
				
					# decide to step forward or back based on value
					if (lnL == float('-inf')):
						population.allele_freq_list[locus_index] = prev_allele_freq
						allele.y = prev_y
						prev_lnL = prev_lnL
						#print("1")
					else:
						value = math.exp(lnL - prev_lnL) * math.exp(prev_y - new_y)
						if (value > 1):
							prev_y = allele.y
							prev_lnL = lnL
							#print("2")
						else:
							random_number = random.random()
							if (random_number < value):
								prev_y = allele.y
								prev_lnL = lnL
								#print("2")
							else:
								population.allele_freq_list[locus_index] = prev_allele_freq
								allele.y = prev_y
								prev_lnL = prev_lnL
								#print("1")
			
			#changes the genotype at a random maternal locus
			for fam in families:
				#print(fam)
				# returns a tuple containing three lists; the first is the list of imputed genotypes 
				# the second the alleles, and the third the allele frequencies at those imputed loci
				imputed = fam.mom.get_imputed_loci(fam.population_name)
				imputed_genotypes = imputed[0]
				imputed_alleles = imputed[1]
				imputed_af_list = imputed[2]
				
				# determines if there are imputed genotypes; if so, new genotypes are tried and
				# likelihood values calculated for each family
				if (len(imputed_genotypes) == 0):
					continue # if no imputed genotypes are present, go on to the next family
				else:
					# one locus at a time is changed in each family; locus chosen randomly
					random_locus = random.randint(0, len(imputed_genotypes) - 1)
					genotype = imputed_genotypes[random_locus]
					allele_list = imputed_alleles[random_locus]
					allele_freq = imputed_af_list[random_locus]
		
					prev_first = genotype.first
					#print(prev_first)
					prev_second = genotype.second
					#print(prev_second)
					prev_fam_lnL = fam.calc_progeny_lnL(population.outcrossing_rate, locus_model)
					#print(prev_fam_lnL)
					# returns a tuple with new maternal alleles = (new_first, new_second)
					new_mom = genotype.impute_new_mom(allele_list, allele_freq, fam.mom.inbreeding_coefficient, random_locus, locus_model)
					new_first = new_mom[0]
					new_second = new_mom[1]
		
					# sets new maternal alleles and calculates family lnL
					genotype.first = new_first
					genotype.second = new_second
					#print(genotype.first)
					#print(genotype.second)
					new_fam_lnL = fam.calc_progeny_lnL(population.outcrossing_rate, locus_model)
					#print(new_fam_lnL)
					if (new_fam_lnL == float('-inf')):
						genotype.first = prev_first
						genotype.second = prev_second
						#print("1")
					else:
						fam_lnL_ratio = (new_fam_lnL - prev_fam_lnL)
						if (fam_lnL_ratio > 0):
							prev_first = genotype.first
							prev_second = genotype.second
							#print("2")
						else:
							rand = random.random()
							value = math.exp(fam_lnL_ratio)
							if (rand < value):
								prev_first = genotype.first
								prev_second = genotype.second
								#print("2")
							else:
								genotype.first = prev_first
								genotype.second = prev_second
								#print("1")
								
		
				for n, genotype in enumerate(fam.mom.genotype_list):
					if str(genotype) not in fam.possible_genotypes[n]:
						fam.possible_genotypes[n].append(str(genotype))
		
				if step > burn_in:
 					if (step % 10) == 0:
 						for n, genotype in enumerate(fam.mom.genotype_list):
 							fam.locus_genotypes[n].append(str(fam.mom.genotype_list[n]))

			if step > burn_in:
				if (step % 10) == 0:
					pop_lnL_list.append(population.calc_pop_lnL(locus_model))
					if(writeOutput3):
						borice_output3.write("%.6f" % population.calc_pop_lnL(locus_model) + "\n")

		endtime = time.asctime()
		print("end time was %s" % endtime)
		
		# main code dealing with file output begins here
		borice_output1.write("Posterior distribution of population inbreeding history:\n")

		for ih_value in range(0, 7):
			ih_count = ih_list.count(ih_value)
			ih_total = len(ih_list)
			ih_percent = float(ih_count)/ih_total
			borice_output1.write("IH value =\t%s\tProportion =\t%.4f\n" % (ih_value, ih_percent))

		borice_output1.write("\nPosterior distributions of t and F:\n")
		borice_output1.write("\nt and F values range from 0 to 1.\nThe below bins represent a range of values greater than or equal to the number listed, but less than the next value.\nExcept bin 1.00, which represents only instances where the t or F value equaled 1.00.\n\n")

		t_percent_list = []
		F_percent_list = []
		bin = 0.0
		while bin < 1.01:
			t_bin_list = []
			for n, t in enumerate(t_list):
				if (t >= bin) and (t < (bin + 0.01)):
					t_bin_list.append(t)
				else:
					continue
			t_count = len(t_bin_list)
			t_percent = float(t_count)/len(t_list)
			t_percent_list.append(t_percent)
			F_bin_list = []
			for n, F in enumerate(F_list):
				if (F >= bin) and (F < (bin + 0.01)):
					F_bin_list.append(F)
				else:
					continue
			F_count = len(F_bin_list)
			F_percent = float(F_count)/len(F_list)
			F_percent_list.append(F_percent)
			borice_output1.write("bin =\t%.2f\tt proportion =\t%.4f\tF proportion =\t%.4f\n" % (bin, t_percent, F_percent))
			bin = bin + 0.01

		t_percent_max = 0.0
		for n, percent in enumerate(t_percent_list):
			if percent > t_percent_max:
				t_percent_max = percent
				t_max = n * 0.01
			else:
				continue
		
		t_list.sort()
		t_lower = int(len(t_list) * 0.025)
		t_lower_percentile = t_list[t_lower]
		t_upper = int(len(t_list) * 0.975)
		t_upper_percentile = t_list[t_upper]
		sum_t = math.fsum(t_list)
		total_t = len(t_list)
		mean_t = (sum_t/total_t)
		borice_output1.write("\nMean t = %.2f; t-max = %.2f; 2.5 percentile = %.2f; 97.5 percentile = %.2f\n" % (mean_t, t_max, t_lower_percentile, t_upper_percentile))

		F_percent_max = 0.0
		for n, percent in enumerate(F_percent_list):
			if percent > F_percent_max:
				F_percent_max = percent
				F_max = n * 0.01
			else:
				continue

		F_list.sort()
		F_lower = int(len(F_list) * 0.025)
		F_lower_percentile = F_list[F_lower]
		F_upper = int(len(F_list) * 0.975)
		F_upper_percentile = F_list[F_upper]
		sum_F = math.fsum(F_list)
		total_F = len(F_list)
		mean_F = (sum_F/total_F)
		borice_output1.write("Mean F = %.2f; F-max = %.2f; 2.5 percentile = %.2f; 97.5 percentile = %.2f\n\n" % (mean_F, F_max, F_lower_percentile, F_upper_percentile))
		
		new_likelihoods = []
		sum_ln_likelihoods = math.fsum(pop_lnL_list)
		total_n = len(pop_lnL_list)
		arithmetic_mean = (sum_ln_likelihoods/total_n)
		borice_output1.write("Ave LL = %s\n" % arithmetic_mean)
# 		for n, ln_likelihood in enumerate(pop_lnL_list):
# 			new_likelihood = Decimal(-(ln_likelihood)).exp()
# 			new_likelihoods.append(new_likelihood)
# 		new_sum = sum(new_likelihoods)
# 		reciprocal_n = math.pow(total_n, -1)
# 		new_reciprocal_n = Decimal(reciprocal_n)
# 		reciprocal_harmonic = new_sum * new_reciprocal_n
# 		harmonic_mean = Decimal(1/reciprocal_harmonic)
# 		borice_output1.write("Harmonic mean of likelihood values = %s\n" % harmonic_mean)
		
		if(writeOutput2):
			borice_output2.write("\nPosterior distributions of maternal inbreeding histories:\n")

		families = sorted(families)
		for fam in families:
			if(writeOutput2):
				borice_output2.write("Maternal Individual %s:\n" % fam.name)
			for ih_value in range(0, 7):
				ih_count = fam.inbreeding_history_list.count(ih_value)
				ih_total = len(fam.inbreeding_history_list)
				ih_percent = float(ih_count)/ih_total
				if(writeOutput2):
					borice_output2.write("IH value =\t%s\tProportion =\t%.4f\n" % (ih_value, ih_percent))

		borice_output1.write("\nPosterior distributions of allele frequencies at each locus:\n")
		borice_output1.write("\nAllele frequency values range from 0 to 1.\nThe below bins represent a range of values greater than or equal to the number listed, but less than the next value.\nExcept bin 1.00, which represents only instances where the allele frequency equaled 1.00.\n\n")
		borice_output1.write("Allele 0 is the null allele. All other alleles are observed in the dataset.\n\n")

		bin = 0.0
		for locus_index, locus in enumerate(all_alleles):
			locus_alleles = all_alleles[locus_index]
			borice_output1.write("Locus %s\n" % (locus_index + 1))
			for n, allele in enumerate(locus_alleles):
				null = locus_model[locus_index]
				if null:
 					borice_output1.write("Allele %s\n" % allele.name)
 					af_list = allele.af_list
 					while bin < 1.01:
 						af_bin_list = []
 						for i, af in enumerate(af_list):
 							if (af >= bin) and (af < (bin + 0.01)):
 								af_bin_list.append(af)
 							else:
 								continue
 						af_count = len(af_bin_list)
 						af_percent = float(af_count)/len(af_list)
 						borice_output1.write("bin =\t%.2f\taf proportion =\t%.4f\n" % (bin, af_percent))
 						bin = bin + 0.01
 					bin = 0.0
 				else:
 					if n == 0: # skips allele zero, which is a dummy allele that remains at a frequency of zero
 						continue
 					else:
						borice_output1.write("Allele %s\n" % allele.name)
						af_list = allele.af_list
						while bin < 1.01:
							af_bin_list = []
							for i, af in enumerate(af_list):
								if (af >= bin) and (af < (bin + 0.01)):
									af_bin_list.append(af)
								else:
									continue
							af_count = len(af_bin_list)
							af_percent = float(af_count)/len(af_list)
							borice_output1.write("bin =\t%.2f\taf proportion =\t%.4f\n" % (bin, af_percent))
							bin = bin + 0.01
						bin = 0.0

		if(writeOutput4):
			borice_output4.write("Posterior distributions for each maternal genotype at each locus in each family.\n")

		for fam in families:
			if(writeOutput4):
				borice_output4.write("Maternal Individual %s:\n" % fam.name)
			for n, locus in enumerate(fam.possible_genotypes):
				if(writeOutput4):
					borice_output4.write("Locus %s\n" % (n + 1))
				possible_genotype_list = fam.possible_genotypes[n]
				locus_list = fam.locus_genotypes[n]
				for i, genotype in enumerate(possible_genotype_list):
					genotype_count = locus_list.count(str(genotype))
					genotype_percent = float(genotype_count)/len(locus_list)
					if(writeOutput4):
						borice_output4.write("possible genotype = \t%s\tproportion =\t%.2f\n" % (genotype, genotype_percent))
		
		#Progress complete
		self.current_step += 1
		
		borice_output1.close()
		if(writeOutput2):
			borice_output2.close()
		if(writeOutput3):
			borice_output3.close()
		if(writeOutput4):
			borice_output4.close()

# required to run BORICE from the command line		
if __name__ == '__main__':
	Borice = BORICE()
	dataFileName = sys.argv[1]
	#locus_model = [bool(int(sys.argv[2])), bool(int(sys.argv[3])), bool(int(sys.argv[4])), bool(int(sys.argv[5])), bool(int(sys.argv[6])), bool(int(sys.argv[7])), bool(int(sys.argv[8])), bool(int(sys.argv[9])), bool(int(sys.argv[10])), bool(int(sys.argv[11]))]
	locus_model = [bool(int(sys.argv[2])), bool(int(sys.argv[3])), bool(int(sys.argv[4]))]
	numSteps = int(sys.argv[5])
	numBurnInSteps = int(sys.argv[6])
	outcrossingRateTuningParam = float(sys.argv[7])
	alleleFreqTuningParam = float(sys.argv[8])
	outcrossingRate = float(sys.argv[9])
	writeOutput2, writeOutput3, writeOutput4 = True, True, True
	try:
		writeOutput2 = bool(int(sys.argv[10]))
		writeOutput3 = bool(int(sys.argv[11]))
		writeOutput4 = bool(int(sys.argv[12]))
	except:
		pass
	sys.stderr.write(repr(dataFileName) + '\n	 ')
	sys.stderr.write(repr(locus_model) + '\n	 ')
	sys.stderr.write(repr(numSteps) + '\n	 ')
	sys.stderr.write(repr(numBurnInSteps) + '\n   ')
	sys.stderr.write(repr(outcrossingRateTuningParam) + '\n   ')
	sys.stderr.write(repr(alleleFreqTuningParam) + '\n	  ')
	sys.stderr.write(repr(outcrossingRate) + '\n	')
	sys.stderr.write(repr(writeOutput2) + '\n	 ')
	sys.stderr.write(repr(writeOutput3) + '\n	 ')
	sys.stderr.write(repr(writeOutput4) + '\n	 ')
		
	Borice.main(dataFileName, locus_model, numSteps, numBurnInSteps, outcrossingRateTuningParam, alleleFreqTuningParam, outcrossingRate, writeOutput2, writeOutput3, writeOutput4)