# This python program was built to take a bunch of
# documents as input and determine how much of a target
# document is similar to the other documents by
# using various algorithms
import os
import time



d = 256
kmpcounter=0

# First the KMP algorith: reference: https://www.geeksforgeeks.org/kmp-algorithm-for-pattern-searching/
# This algorithm is for checking single character matching and will be weighed the least
def KMPCheck(pat, txt):
    M = len(pat)
    N = len(txt)
    kmpcounter = 0
    # longest prefix suffix:
    lps = [0] * M
    j = 0  # pattern index

    computeLPSArray(pat, M, lps)

    i = 0;  # txt index
    while i < N:
        if pat[j] == txt[i]:
            i += 1
            j += 1
        if j == M:
            kmpcounter+=1
            print("Found pattern at index " + str(i - j))
            j = lps[j - 1]

        elif i < N and pat[j] != txt[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return kmpcounter

def computeLPSArray(pat, M, lps):
    len = 0  # length of previous lps
    lps[0]
    i = 1

    while i < M:
        if pat[i] == pat[len]:
            len += 1
            lps[i] = len
            i += 1
        else:
            if len != 0:
                len = lps[len - 1]

            else:
                lps[i] = 0
                i += 1


# Run the LCSS algorithm for each paragraph from the test file and existing corpora.
# Given two sequences, find the length of longest subsequence present in both of them.
# A subsequence is a sequence that appears in the same relative order, but not necessarily
# contiguous. For example, “abc”, “abg”, “bdf”, “aeg”, ‘”acefg”, .. etc are subsequences of “abcdefg”.
# https://www.geeksforgeeks.org/longest-common-subsequence-dp-4/





def lcs(X, Y):
    # find the length of the strings
    m = len(X)
    n = len(Y)

    # declaring the array for storing the dp values
    L = [[None] * (n + 1) for i in range(m + 1)]

    """Following steps build L[m+1][n+1] in bottom up fashion 
    Note: L[i][j] contains length of LCS of X[0..i-1] 
    and Y[0..j-1]"""
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 or j == 0:
                L[i][j] = 0
            elif X[i - 1] == Y[j - 1]:
                L[i][j] = L[i - 1][j - 1] + 1
            else:
                L[i][j] = max(L[i - 1][j], L[i][j - 1])

                # L[m][n] contains the length of LCS of X[0..n-1] & Y[0..m-1]
    return L[m][n]



# end of function lcs
# Now the Rabin-Karp Algorithm, which checks for entire words that match:
# https://www.geeksforgeeks.org/rabin-karp-algorithm-for-pattern-searching/

def RBKCheck(pat, txt, q):  # q must be a prime number
    rbkcounter = 0
    M = len(pat)
    N = len(txt)
    i = 0
    j = 0
    p = 0
    t = 0  # pattern hash
    h = 1  # txt hash
    for i in range(M - 1):
        h = (h * d) % q
    for i in range(M):
        p = (d * p + ord(pat[i])) % q
        t = (d * t + ord(txt[i])) % q

    for i in range(N - M + 1):
        if p == t:
            for j in range(M):
                if txt[i + j] != pat[j]:
                    break
            j += 1
            if j == M:
                print("pattern found at index  "+ str(i))
                rbkcounter+=1
        if i < N - M:
            t = (d * (t - ord(txt[i]) * h) + ord(txt[i + M])) % q

            if t < 0:
                t = t + q
    return rbkcounter

# here is where the main program will run

print("here is the current directory: " + os.getcwd())

#
#print('Then cost is ',time_end-time_start,'s' )

#txt = "AABAACAADAABAABA"

#at="AABA"
#txt = "GEEKS FOR GEEKS"
#pat = "GEEK"
q = 101
#RBKCheck(pat, txt, q)
#print(RBKCheck(pat, txt, q))
#X="AGGTAB"
#Y="GXTXAYB"


#print(KMPCheck(pat, txt))
#KMPCheck(pat, txt)+1
#print(KMPCheck(pat, txt))



#first example of input
time_start=time.time()
writtenfile = open("example.txt", "r")
newfile = open("newfile.txt", "r")

txt= writtenfile.read()  # this text file now holds the entire string from example.txt
ntxt = newfile.read()
lcslength=lcs(txt, ntxt)

similarity=(KMPCheck(ntxt, txt)*0.3)+(lcslength/len(txt)*1)+(RBKCheck(ntxt, txt, q)*0.10)
print("The similarity of newfile is", similarity)
if similarity>1:
    print("It is plagiarized:")
elif similarity>0.5:
    print("It is alarmingly similar.")
else:
    print("It is likely just referenced.")
time_end=time.time()
print('Totally cost of inputsize',len(ntxt),'is',time_end-time_start)
print('\n')

#second example of input
time_start=time.time()
w1= open("e1.txt", "r")
n1 = open("n1.txt", "r")

txt1= w1.read()  # this text file now holds the entire string from example.txt
ntxt1 = n1.read()

lcslength=lcs(txt1, ntxt1)

similarity=(KMPCheck(ntxt1, txt1)*0.3)+(lcslength/len(txt1)*1)+(RBKCheck(ntxt1, txt1, q)*0.10)
print("The similarity is of second example", similarity)
if similarity>1:
    print("It is plagiarized:")
elif similarity>0.5:
    print("It is alarmingly similar.")
else:
    print("It is likely just referenced.")
time_end=time.time()
print('Totally cost of inputsize',len(ntxt1),'is',time_end-time_start)
print('\n')


#third example of input
time_start=time.time()
w2= open("e2.txt", "r")
n2 = open("n2.txt", "r")

txt2= w2.read()  # this text file now holds the entire string from example.txt
ntxt2 = n2.read()

lcslength=lcs(txt2, ntxt2)

similarity=(KMPCheck(ntxt2, txt2)*0.3)+(lcslength/len(txt2)*1)+(RBKCheck(ntxt2, txt2, q)*0.10)
print("The similarity is of third example", similarity)
if similarity>1:
    print("It is plagiarized:")
elif similarity>0.5:
    print("It is alarmingly similar.")
else:
    print("It is likely just referenced.")
time_end=time.time()
print('Totally cost of inputsize',len(ntxt2),'is',time_end-time_start)