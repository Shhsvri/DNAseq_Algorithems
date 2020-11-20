# Boyer Moore implementation

from bm_preproc import *

p = 'TCAA'
p_bm = BoyerMoore(p)
# print(p_bm.bad_character_rule(2, 'T'))

p = 'ACTA'
p_bm = BoyerMoore(p)
p_bm.good_suffix_rule(0)

p = 'ACAC'
p_bm = BoyerMoore(p)
p_bm.match_skip()

def boyer_moore(p, p_bm, t):
    i = 0
    occurrences = []
    while i < len(t)-len(p)+1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if not p[j] == t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences


t= 'GCTACGATCTAGAATCTA'
p= 'TCTA'
p_bm = BoyerMoore(p)
# print(boyer_moore(p,p_bm, t))
