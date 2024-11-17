import sys
import random
from collections import defaultdict, deque

def kmer_count(sequence, randomized, k):
  print(f"{'K-mer':<10}{'Original':<10}{'Randomized':<10}")
  print("-" * 30)
  seq_kmers = defaultdict(int)
  for i in range(len(sequence) - k + 1):
    kmer = sequence[i:i+k]
    seq_kmers[kmer] += 1

  rand_kmers = defaultdict(int)
  for i in range(len(randomized)-k+1):
    kmer = randomized[i:i+k]
    rand_kmers[kmer] += 1

  all_kmers = set(seq_kmers.keys()).union(rand_kmers.keys())
  for kmer in all_kmers:
    print(f"{kmer:<10}{seq_kmers[kmer]:<10}{rand_kmers[kmer]:<10}")

def random_euler(filename, k):
  try:
    with open(filename, 'r') as file:
      header = file.readline().strip()
      contents = [line.strip() for line in file.readlines()]
      sequence = ''.join(contents)
  except:
    print("The file does not exist")
    return

  kmers = []
  edges = []
  graph = defaultdict(list)
  edge_count = defaultdict(int)
  start = sequence[0:k]
  ending = sequence[-k:]
  total_edges = 0
  for i in range(len(sequence) - k + 1):
    kmer = sequence[i:i+k]
    prefix = kmer[:-1]
    suffix = kmer[1:]
    edges.append((prefix, suffix))
    #print(f"{prefix}, {suffix}")
    edge_count[(prefix, suffix)] += 1
    graph[prefix].append(suffix)
    kmers.append(kmer)
  og_edge_count = edge_count # in case restarting is necessary
  nodes = list(graph.keys())
  curr = start[1:]
  randomized = start
  used = []
  backtraced = []
  total_edges = sum(len(edges) for edges in graph.values())
  fits = False
  while not fits :
    while len(randomized) < len(sequence) - 1:
      out = graph[curr]
      out1 = [nex for nex in out if edge_count[(curr, nex)] > 0 and nex != ending[1:]] #checks that it's not the last kmer and that it's not already used
      out = out1
      if out:
        nex = random.sample(out, 1)[0]
        edge = (curr, nex)
        edge_count[edge] -= 1
        used.append(edge)
        randomized += nex[-1]
        curr = nex
      else:
        if len(used) == 0: #if backtracing goes too far, start sequence over. theoretically shouldn't happen
          curr = start[1:]
          backtraced = []
          used = []
          randomized = start
          edge_count = og_edge_count
        else:
          last = used.pop()
          randomized = randomized[:-1]
          edge_count[last] += 1
          while last in backtraced and len(used) > 0: #if the edge has already been backtraced, go even further
            last = used.pop()
            edge_count[last] += 1
            randomized = randomized[:-1]
          curr = last[0]
          backtraced.append(last)
          if len(used) == 0:
            backtraced = []

    lasts = graph[curr]
    if lasts:
      if ending[1:] in lasts:
        fits = True
      else:
        if len(used) == 0: #seemed to occur when there is a cycle earlier in the implementation process
          curr = start[1:]
          backtraced = []
          randomized = start
          edge_count = og_edge_count
        else:
          lastone = used.pop()
          randomized = randomized[:-1]
          while lastone in backtraced and len(used) > 0:
            lastone = used.pop()
            randomized = randomized[:-1]
          curr = lastone[0]
          backtraced.append(lastone)
    else: #if lasts has nothing, theoretically should never go here
      used = []
      curr = start[1:]
      backtraced = []
      randomized = start
      edge_count = og_edge_count
  randomized += ending[-1]
  output_filename = "randomized_seq.fasta"
  output_header = header + "_randomized"

  with open(output_filename, "w") as out:
    out.write(output_header + "\n")
    out.write(randomized)

  kmer_count(sequence, randomized, k)



def randomize(filename, k):
  try:
    with open(filename, 'r') as file:
      header = file.readline().strip()
      contents = [line.strip() for line in file.readlines()]
      sequence = ''.join(contents)
  except:
    print("The file does not exist")
    return

  kmers = []

  for i in range(len(sequence) - k + 1):
    kmer = sequence[i:i+k]
    kmers.append(kmer)

  total_kmers = len(kmers)

  randomized = random.sample(kmers, 1)[0]
  kmers.remove(randomized)

  possible = []
  subk = k - 1
  used = []
  used.append(randomized)

  while len(kmers) != 0:
    for kmer in kmers:
      if kmer[0:subk] == randomized[-subk:]:
        possible.append(kmer)
    if len(possible) == 0: #if there are no kmers with the ending substring
      possible = kmers
      chosen = random.sample(possible, 1)[0] #looks for a random one and adds the last letter to the sequence --> different kmer content
      randomized += chosen[subk];
      used.append(chosen)
      kmers.remove(chosen)
      possible = []
    else:
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

