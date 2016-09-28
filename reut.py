from suffix import suffix_tree, gst


all_res = [['A.509', 'A.512', 'A.512', 'A.514', 'A.515', 'A.516', 'A.517', 'A.523', 'A.524', 'A.527', 'A.527', 'A.528', 'A.121', 'A.124', 'A.125', 'A.126', 'A.133', 'A.134', 'A.136', 'A.136', 'A.137', 'A.138', 'A.221', 'A.228', 'A.242', 'A.265', 'A.497', 'A.500', 'A.512', 'A.512', 'A.514', 'A.516', 'A.517', 'A.523', 'A.524', 'A.527', 'A.527', 'A.528', 'A.531', 'A.497', 'A.504', 'A.510', 'A.511', 'A.512', 'A.515', 'A.516', 'A.517', 'A.517', 'A.519', 'A.524', 'A.527', 'A.528', 'A.531', 'A.107', 'A.107', 'A.111', 'A.112', 'A.115', 'A.119', 'A.119', 'A.119', 'A.122', 'A.128', 'A.132', 'A.133', 'A.242', 'A.249', 'A.454', 'A.459', 'A.459', 'A.460', 'A.463', 'A.464', 'A.466', 'A.466', 'A.526', 'A.527', 'A.527', 'A.530', 'A.530', 'A.534', 'A.118', 'A.119', 'A.119', 'A.120', 'A.121', 'A.121', 'A.122', 'A.138', 'A.141', 'A.142', 'A.270', 'A.270', 'A.125', 'A.126', 'A.126', 'A.130', 'A.133', 'A.133', 'A.134', 'A.136', 'A.138', 'A.223', 'A.227', 'A.111', 'A.120', 'A.121', 'A.122', 'A.122', 'A.124', 'A.125', 'A.125', 'A.128', 'A.140', 'A.143', 'A.146', 'A.217', 'A.218', 'A.219', 'A.221', 'A.222', 'A.222', 'A.246', 'A.252', 'A.281', 'A.282'], ['A.119', 'A.121', 'A.122', 'A.142', 'A.142', 'A.143', 'A.143', 'A.143', 'A.146', 'A.148', 'A.149', 'A.102', 'A.105', 'A.107', 'A.283', 'A.296', 'A.296', 'A.298', 'A.298', 'A.299', 'A.300', 'A.121', 'A.124', 'A.141', 'A.142', 'A.143', 'A.143', 'A.146', 'A.147', 'A.149', 'A.149', 'A.151', 'A.151', 'A.121', 'A.122', 'A.126', 'A.126', 'A.142', 'A.143', 'A.143', 'A.143', 'A.146', 'A.147', 'A.149', 'A.149', 'A.150', 'A.151', 'A.307', 'A.309', 'A.318', 'A.318', 'A.319', 'A.321', 'A.324', 'A.325', 'A.325', 'A.136', 'A.137', 'A.140', 'A.141', 'A.144', 'A.147', 'A.148', 'A.150', 'A.151', 'A.153', 'A.213', 'A.217', 'A.217', 'A.221', 'A.222', 'A.224', 'A.225', 'A.226', 'A.251', 'A.254', 'A.211', 'A.217', 'A.218', 'A.221', 'A.221', 'A.250', 'A.251', 'A.252', 'A.205', 'A.217', 'A.218', 'A.218', 'A.219', 'A.228', 'A.230', 'A.233', 'A.236', 'A.239', 'A.345', 'A.349', 'A.350', 'A.352', 'A.361', 'A.370', 'A.370', 'A.370', 'A.380', 'A.511', 'A.511']]
data = [" ".join(all_res[0]), " ".join(all_res[1])]

st = gst.STree(all_res)
print "!!"
a,b=st.lcs() # "abc"


tree0 = suffix_tree.SuffixTree(all_res[0])
# tree1 = suffix_tree.SuffixTree(all_res[1])
index_of_need = tree0.find_substring(['A.121', 'A.124'])
# print index_of_need
# def long_substr(data):
#     substr = ''
#     if len(data) > 1 and len(data[0]) > 0:
#         for i in range(len(data[0])):
#             for j in range(len(data[0])-i+1):
#                 if j > len(substr) and is_substr(data[0][i:i+j], data):
#                     substr = data[0][i:i+j]
#     return substr

# def is_substr(find, data):
#     if len(data) < 1 and len(find) < 1:
#         return False
#     for i in range(len(data)):
#         if find not in data[i]:
#             return False
#     return True


# print long_substr(['Oh, hello, my friend.',
#                    'I prefer Jelly Belly beans.',
#                    'When hell freezes over!'])

def issublist(a, b):
	# import pdb
	# pdb.set_trace()
	return " ".join(a) in " ".join(b)
	
def long_substr(data):
    max_len = 0
    substr_by_len = dict()
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                # if j > len(substr) and all(issublist(data[0][i:i+j], x) for x in data):
                if j > 0 and all(issublist(data[0][i:i+j], x) for x in data):
                    substr = " ".join(data[0][i:i+j])
                    if j not in substr_by_len:
                    	substr_by_len[j] = set([substr])
                    else:
                    	substr_by_len[j].add(substr)
                    if j > max_len:
                    	max_len = j

                    # pdb.set_trace()
                    

    return max_len, substr_by_len

max_len, substr_by_len= long_substr(all_res )
# print max_len
print len(a)
print len(substr_by_len[max_len])
print a
print b
print "____"
print substr_by_len[max_len]
print a==b
# long_substr(all_res)
