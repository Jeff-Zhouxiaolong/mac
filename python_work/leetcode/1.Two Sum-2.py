#改代码在第6行“for j , b in enumerate(numbers[i + 1 - len(numbers)]):”
#其中enumerate()里面的内容觉得有问题，但是该方法解这道题完全没问题
class Solution:
    def twoSum(self, numbers, target):
        for i, a in enumerate(numbers):
            for j , b in enumerate(numbers[i + 1 - len(numbers)]):
                if a + b == target:
                    return [i, j + i + 1]
        return [-1, -1]


