import bm_preproc as bm
from boyer_moore import boyer_moore

def approximate_match(p, t, n):
    segment_length = round(len(p) / (n+1))
    all_matches = set()
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        p_bm = bm.BoyerMoore(p[start:end], alphabet='ACGT')
        matches = boyer_moore(p[start:end], p_bm, t)
        
        for m in matches:
            if m < start or m-start-len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0,start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches >n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches)

p = 'AACTTG'
t = 'CACTTAATTTG'

print(approximate_match(p, t, 2))
