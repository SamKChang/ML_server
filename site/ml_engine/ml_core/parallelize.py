#!/usr/bin/env python
# python parallel wrapper for simply defined functions
# it returns a list of output corresponding to each input entry

import multiprocessing as mp
import operator
from compiler.ast import flatten

def parallelize(target_function,
		input_list,
		threads,
		**kwargs):
	"""
	target_function is implemented in a general way
	supposely any function would work
	But it could break down if target_function assumes some convoluted data structure

	input_list is a list of list. 
	Each input entry should be wrapped properly as a list 
	**kwargs can be passed py passing dictionary

	Example:
	# a toy target function
	def f(a, b, **kwargs):
	if 'factor' in kwargs:
	factor = kwargs['factor']
	else:
	factor = 1
	return a + b*factor

	input_list = [[i,j,{'factor':3}] for i in range(10) for j in range(10)]

	out_list = parallelize(f, input_list, 2, block_size=10)
	"""

	if 'block_size' in kwargs:
		block_size = kwargs['block_size']
	else:
	    # default block_size: each thread gets 3 queues
		block_size = len(input_list)/(threads*3)
    
  #############################################
  # runing target function of a single thread #
  #############################################
	def run_jobs(q_in, q_out):
		for inps in iter(q_in.get, None):
			ind = inps[-1]    # index of job
			inps = inps[:-1]  # actual input sequence
			out = []
			for args in inps:
				if type(args[-1]) == dict: # check known args input
					kwargs = args[-1]
					args = args[:-1]
					out.append(target_function(*args, **kwargs))
				else:
					out.append(target_function(*args))
			if out != None:
				q_out.put([out, ind]) # output result with index
			q_in.task_done()
		q_in.task_done() # task done for 'None' if q_in finished
	  ###### end of single thread definition ######
  
  # devide input_list into chunks according to block_size
	def chunks(_list, _size):
		for i in xrange(0, len(_list), _size):
			yield _list[i:i+_size]
	input_block = list(chunks(input_list, block_size))

  # setup empty queue
	output_stack = []
	output = []
	qinp = mp.JoinableQueue()
	qout = mp.Queue()

  # start process with empty queue
	for thread in range(threads):
		p =  mp.Process(target=run_jobs, args=(qinp, qout))
		p.daemon = True # necessary for terminating finished thread
		p.start()

    # put I/O data into queue for parallel processing
	index = range(len(input_block))
	for ind, inps in zip(index, input_block):
		inps.append(ind) # append inp index
		qinp.put(inps)   # put inp to input queue
	qinp.join()       # wait for jobs to finish

  # 'while not queue.empty' is NOT reliable
	if not qout.empty():
		for i in range(len(input_block)):
			output_stack.append(qout.get())

	if len(output_stack)>0:
		# sort/restructure output according to input order
		output_stack = sorted(output_stack, key=operator.itemgetter(1))
		# loop though all input for corresponding output
		for data_out in output_stack:
			# if output is list of class, in-line iteration doesn't work
			output.append(data_out[0])
		return flatten(output)
