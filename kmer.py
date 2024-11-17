import sys
import random
from collections import defaultdict, deque

def kmer_count(sequence, randomized, k):
  print(f"{'K-mer':<10}{'Original':<10}{'Randomized':<10}")
  print("-" * 30)
  seq_kmers = defaultdict(int)
  for i in range(len(sequence) - k + 1):
    kmer = sequence[i:i+k]
   # if kmer in kmers:
    #  kmers[kmer] += 1
    #else:
    #  kmers[kmer] = 1
    seq_kmers[kmer] += 1

  rand_kmers = defaultdict(int)
  for i in range(len(randomized)-k+1):
    kmer = randomized[i:i+k]
    rand_kmers[kmer] += 1

  all_kmers = set(seq_kmers.keys()).union(rand_kmers.keys())
  for kmer in all_kmers:
    print(f"{kmer:<10}{seq_kmers[kmer]:<10}{rand_kmers[kmer]:<10}")

def random_euler(filename, k):
  #only works for k > 1
  try:
    with open(filename, 'r') as file:
      header = file.readline().strip()
      contents = [line.strip() for line in file.readlines()]
      sequence = ''.join(contents)
      #print(sequence)
  except:
    print("The file does not exist")
    return

 # kmers = {}
  kmers = []
  edges = []
  graph = defaultdict(list)
  edge_count = defaultdict(int)
  #print("Original Kmers")
  #kmer_count(sequence, k)
  start = sequence[0:k]
  #print(f"start:{start}")
  ending = sequence[-k:]
  total_edges = 0
  #sequence = sequence[1:-1] #Starting and ending are never added to the graph
  for i in range(len(sequence) - k + 1):
    kmer = sequence[i:i+k]
   # if kmer in kmers:
    #  kmers[kmer] += 1
    #else:
    #  kmers[kmer] = 1
    prefix = kmer[:-1]
    suffix = kmer[1:]
    edges.append((prefix, suffix))
    #print(f"{prefix}, {suffix}")
    edge_count[(prefix, suffix)] += 1
    graph[prefix].append(suffix)
    kmers.append(kmer)
  '''
  for i in range(len(edges)):
    for j in range(len(edges)):
      if edges[i] == edges[j] and i != j:
        print(f"Multiple:{edges[i]}")
  '''
  #return
  '''
  for kmer in kmers:
    suffix = kmer[1:]
    for k in kmer:
      if k != kmer:
        prefix = k[:-1]
        if prefix == suffix:
          graph[kmer].append(k)
          total_edges += 1
  '''
  
  nodes = list(graph.keys())
  #print(nodes)
  '''
  for node in nodes:
    if node == ending[1:]:
      print("same")
    for suffix in graph[node]:
      if suffix == ending[1:]:
        print("suffix found")
  '''

  #return
  #start = random.sample(nodes, 1)[0]
  curr = start[1:]
  randomized = start
  used = []
  #used.append(start)
  #used.append(ending)
  backtraced = []
  total_edges = sum(len(edges) for edges in graph.values())
  #print(total_edges)
  iteration_cap = 20
  backtracks = 0
  fits = False
  multiple = curr
  while not fits :
    while len(randomized) < len(sequence) - 1:
      out = graph[curr]
      #print(out)
      #if (len(out) == 0) :
      #  print("out")

      out1 = [nex for nex in out if edge_count[(curr, nex)] > 0 and nex != ending[1:]]
      out = out1
      if out:
        #print("found edges")
        if len(out) > 1:
          multiple = curr
          #print(f"multiple:{multiple}")
        nex = random.sample(out, 1)[0]
        #print(nex)
        edge = (curr, nex)
        #print(len(out))
        edge_count[edge] -= 1
        #if edge_count[edge] == 0:
        used.append(edge)
        randomized += nex[-1]
        curr = nex
        #print(f"curr:{curr}")
        #print(f"randomized:{randomized}")
      else:
        if len(used) == 0:
          #print("too many backtracks, starting over")
          curr = start[1:]
          backtraced = []
          used = []
          randomized = start
          backtracks = 0
          #return
        else:
          #print("retrace")
          #print(f"before: {randomized}")
          #print(used)
          #print(out)
          last = used.pop()
          randomized = randomized[:-1]
          edge_count[last] += 1
          #for i in range(backtracks):
          
          while last in backtraced and len(used) > 0: #if the edge has already been traced, go even further
            #curr = last[0] #goes back to last element
            #backtraced.remove(last)
            last = used.pop()
            edge_count[last] += 1
            randomized = randomized[:-1]
            #backtracks += 1
            #last = used.pop()
          '''
          while last[0] != multiple: #backtracks to last root with multiple children
            last = used.pop()
            randomized = randomized[:-1]
            edge_count[last] += 1
          '''
          curr = last[0]
          backtraced.append(last)
          backtracks += 1
          if len(used) == 0:
            backtraced = []
          '''
          curr = start[1:]
          backtraced = []
          used = []
          randomized = start
          backtracks = 0
          '''
          #print(f"after: {randomized}")
          #print(f"curr: {curr}")
          #print(f"multiple, retrace: {multiple}")
      #print(randomized)

    #print(randomized)
    #print("out of inner while")
    lasts = graph[curr]
    if lasts:
      if ending[1:] in lasts:
        fits = True
      else:
        if len(used) == 0: #seems to occur when there is a cycle
          curr = start[1:]
          backtraced = []
          randomized = start
          #print("cycle detected")
          #print(used)
        else:
          lastone = used.pop()
          randomized = randomized[:-1]
          while lastone in backtraced and len(used) > 0:
            #print("lastone loop")
            lastone = used.pop()
            randomized = randomized[:-1]
            #lastone = used.pop()
          curr = lastone[0]
          backtraced.append(lastone)
          #print("have to repeat")
          #print(f"randomized:{randomized}")
          #print(f"used:{used}")
          #print(f"backtraced:{backtraced}")
          #print(f"lastone:{lastone}")
    else: #if lasts has nothing
      #print("lasts has nothing, starting over")
      used = []
      curr = start[1:]
      backtraced = []
      randomized = start
  randomized += ending[-1]
  #print(ending)
  output_filename = "randomized_seq.fasta"
  output_header = header + "_randomized"

  with open(output_filename, "w") as out:
    out.write(output_header + "\n")
    out.write(randomized)

  #print("Final kmers")
  kmer_count(sequence, randomized, k)



def randomize(filename, k):
  try:
    with open(filename, 'r') as file:
      header = file.readline().strip()
      contents = [line.strip() for line in file.readlines()]
      sequence = ''.join(contents)
      #print(sequence)
  except:
    print("The file does not exist")
    return

 # kmers = {}
  kmers = []

  for i in range(len(sequence) - k + 1):
    kmer = sequence[i:i+k]
   # if kmer in kmers:
    #  kmers[kmer] += 1
    #else:
    #  kmers[kmer] = 1
    kmers.append(kmer)
  
  #print(kmers)
  
  total_kmers = len(kmers)

  randomized = random.sample(kmers, 1)[0]
  kmers.remove(randomized)
  #used.append(randomized)

  #used = []
  possible = []
  subk = k - 1
  used = []
  used.append(randomized)
  iterations = 0
  prev_unique = 1
  removes = 1

  while len(kmers) != 0:
    for kmer in kmers:
      if kmer[0:subk] == randomized[-subk:]:
        possible.append(kmer)
    if len(possible) == 0: #if there are no kmers with the ending substring
      #print(randomized)
      #print(kmers)
      #print(f"Nothing in possible {randomized[-subk:]}")
      possible = kmers
      chosen = random.sample(possible, 1)[0]
      randomized += chosen[subk];
      used.append(chosen)
      kmers.remove(chosen)
      possible = []
    else:
      #print("possible")
      chosen = random.sample(possible, 1)[0]
      randomized += chosen[subk];
      used.append(chosen)
      kmers.remove(chosen)
      possible = []

  output_filename = "randomized_seq.fasta"
  output_header = header + "_randomized"

  with open(output_filename, "w") as out:
    out.write(output_header + "\n")
    out.write(randomized)

  kmer_count(sequence, randomized, k)

def main():
  if len(sys.argv) != 4:
    print("Usage: python kmer.py <file_name>.fasta <k> <s/e>")
    print("If using Euler: DNA sequences only work with k > 5 and protein sequences only work with k > 2")
    print("Otherwise, program will be stuck in an infinite loop")
  elif sys.argv[3] == "s":
    randomize(sys.argv[1], int(sys.argv[2]))
  elif sys.argv[3] == "e":
    random_euler(sys.argv[1], int(sys.argv[2]))
  else:
    print("Usage: python kmer.py <file_name>.fasta <k> <s/e>")
    print("If using Euler: DNA sequences only work with k > 5 and protein sequences only work with k > 2")
    print("Otherwise, program will be stuck in an infinite loop")


if __name__ == "__main__":
  main()

