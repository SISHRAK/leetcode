from itertools import *
class Solution:
    def sortByBits(self, arr: List[int]) -> List[int]:
        def F(num):
            bit_count = bin(num).count('1')
            return (bit_count, num)
        arr.sort(key = F)
        return arr  
    def countNumbersWithUniqueDigits(self, n: int) -> int:
        ans = 1
        temp = 1
        if(n == 1):
            return 10
        if(n == 0):
            return 1
        for i in range (1, n+1):
            ans = 9*temp + ans
            temp = temp * (10 - i)
        return ans   
