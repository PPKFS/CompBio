import matplotlib.pyplot as plt
import math

window_size = 10000
sequence = []
def get_next_sequence():
	for i in range(len(sequence)-(window_size+1)):
		yield sequence[i:i+window_size]

V = []

def viterbi(obs, states, start, trans, emit):
	#list (time-indexed) of dictionaries (state->prob)
	V = [{}]

	#not sure why there's 2
	newPath = {}
	path = {}
	#seed it with initial probabilities
	initial_probs = [math.log(emit[s][obs[0]]*start[s], 2) for s in states]
	V[0]['h'] = initial_probs[0]
	V[0]['l'] = initial_probs[1]
	#path to start is empty
	path['h'] = []
	path['l'] = []

	for i in range(1, len(obs)):
		#add an empty dict for t=i
		V.append({})
		for y in states:
			#for every state, find the maximum probability of reaching this state from every other state at t=i-1, and the state that gives you that max probability
			(prob, state) = max((math.log(emit[y][obs[i]], 2) + math.log(trans[s][y], 2) + V[i-1][s], s) for s in states)
			#at time i, the probability of being in y is prob
			V[i][y] = prob
			#debug
			print(y, prob)
			#the new path to y is the path to the previous best state, plus y
			newPath[y] = path[state] + [y]
		path = newPath

	#the sequence probability is at t=end, find the highest chance and the corresponding state
	(prob, state) = max((V[len(obs)-1][s], s) for s in states)
	#path[state] is then the previous best state for state, and then the previous best state for the previous best state, recurse
	return (path[state], prob)




if __name__ == "__main__":
	(path, prob) = viterbi('GGATCGC', ['h', 'l'], {'h':0.5, 'l':0.5}, 
		{'h':{'h':0.5, 'l':0.5}, 'l':{'h':0.45, 'l':0.55}},
		{'h': {'T':0.23, 'C':0.29, 'A':0.21, 'G':0.27}, 
		'l': {'T':0.3, 'C':0.21, 'A':0.29, 'G':0.20}})
	print(path, prob)
	with open("phage.txt", "r") as f:
		f.readline()
		for line in f:
			for letter in line:
				if letter != '\n':
					sequence.append(letter)
	dinucs = [(sum(1 if (window[i] == 'C' or window[i] == 'G') else 0 for i in range(len(window)-1)) / (len(window)-1)) for window in get_next_sequence()]
	(path, _) = viterbi(sequence, ['h', 'l'], {'h':0.5, 'l':0.5}, 
		{'h':{'h':0.9997, 'l':0.0003}, 'l':{'h':0.0002, 'l':0.9998}},
		{'h': {'T':0.2079, 'C':0.2478, 'A':0.2459, 'G':0.2984}, 
		'l': {'T':0.3237, 'C':0.2080, 'A':0.2698, 'G':0.1985}})
	#for window in get_next_sequence():	
	#	expected = sum(1 if (x == 'C') else 0 for x in window) * sum(1 if (x == 'G') else 0 for x in window) / len(window)
	#	observed = sum(1 if window[i:i+2] == ['C', 'G'] else 0 for i in range(len(window)-1))
	#	dinucs.append(observed/expected)*/
	#print(dinucs)
	plt.figure(1)
	plt.title('Level of CG content')
	plt.plot(range(len(dinucs)), [0.5 for _ in range(len(dinucs))])
	plt.plot(range(len(dinucs)), dinucs)
	isHigh = False
	start = 0
	for i in range(len(path)):
		if path[i] == 'h':
			if not isHigh:
				isHigh = True
				start = i
		else:
			if isHigh:
				plt.axvspan(start, i, color='red', alpha=0.5)
				print(start, i)
				isHigh = False
	if isHigh:
		plt.axvspan(start, len(path)-1,  color='red', alpha=0.5)
	plt.show()