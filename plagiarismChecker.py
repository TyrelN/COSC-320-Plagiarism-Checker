# This python program was built to take a bunch of 
# documents as input and determine how much of a target 
# document is similar to the other documents by 
# using various algorithms
import os
d = 256


#First the KMP algorith: reference: https://www.geeksforgeeks.org/kmp-algorithm-for-pattern-searching/
#This algorithm is for checking single character matching and will be weighed the least
def KMPCheck(pat, txt):
    M = len(pat)
    N = len(txt)

    #longest prefix suffix:
    lps = [0]*M
    j = 0 #pattern index

    computeLPSArray(pat, M, lps)

    i=0; #txt index
    while i < N:
        if pat[j] == txt[i]:
            i += 1
            j += 1
        if j == M:
            print("Found pattern at index" + str(i-j))
        elif i < N and pat[j] != txt[i]:
            if j != 0:
                j= lps[j-1]
            else:
                i+=1

def computeLPSArray(pat, M, lps):
    len = 0 #length of previous lps
    lps[0] 
    i = 1

    while i < M:
        if pat[i] == pat[len]:
            len += 1
            lps[i] = len
            i+= 1
        else:
            if len!= 0:
                len = lps[len-1]

            else:
                lps[i] = 0
                i += 1

#Run the LCSS algorithm for each paragraph from the test file and existing corpora.
#Given two sequences, find the length of longest subsequence present in both of them.
#A subsequence is a sequence that appears in the same relative order, but not necessarily 
#contiguous. For example, “abc”, “abg”, “bdf”, “aeg”, ‘”acefg”, .. etc are subsequences of “abcdefg”. 
#https://www.geeksforgeeks.org/longest-common-subsequence-dp-4/               

def lcs(X, Y, m, n): 

    if m == 0 or n == 0: 
       return 0; 
    elif X[m-1] == Y[n-1]: 
       return 1 + lcs(X, Y, m-1, n-1); 
    else: 
       return max(lcs(X, Y, m, n-1), lcs(X, Y, m-1, n));                 

def lcs(X, Y):
    # find the length of the strings 
    m = len(X)
    n = len(Y)

    # declaring the array for storing the dp values 
    L = [[None] * (n + 1) for i in xrange(m + 1)]

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

#Now the Rabin-Karp Algorithm, which checks for entire words that match:
#https://www.geeksforgeeks.org/rabin-karp-algorithm-for-pattern-searching/

def RBKCheck(pat, txt, q):# q must be a prime number
    M = len(pat)
    N = len(txt)
    i = 0
    j = 0
    p = 0
    t = 0 #pattern hash
    h = 0 #txt hash 
    for i in range(M-1):
        h = (h*d)%q
    for i in range(M):
        p = (d*p + ord(pat[i]))%q
        t = (d*t + ord(txt[i]))%q

    for i in range(N-M+1):
        if p == t:
            for j in range(M):
                if txt[i + j] != pat[j]:
                    break
            j+=1
            if j==M:
                print("pattern found at index " + str(i))
        
        if i < N-M:
            t = (d*(t-ord(txt[i])*h) + ord(txt[i + M]))%q

            if t<0:
                t = t + q


#here is where the main program will run
print("here is the current directory: " + os.getcwd())

testedfile = open("example.txt", "r")
testfile = open("test.txt", "r")

text1 = testedfile.read()  # this text file now holds the entire string from example.txt
text2 = testedfile.read()
print(text1)

time_start=time.time()

time_end=time.end()
print('Then cost is ',time_end-time_start,'s' )
