'''
给定一个字符串，找出不含有重复字符的最长子串的长度。

示例 1:

输入: "abcabcbb"
输出: 3 
解释: 无重复字符的最长子串是 "abc"，其长度为 3。
示例 2:

输入: "bbbbb"
输出: 1
解释: 无重复字符的最长子串是 "b"，其长度为 1。
示例 3:

输入: "pwwkew"
输出: 3
解释: 无重复字符的最长子串是 "wke"，其长度为 3。
     请注意，答案必须是一个子串，"pwke" 是一个子序列 而不是子串。
'''
class Solution:
    def lengthOfLongestSubstring(self, s):
        unique_chars = set([])
        j = 0
        n = len(s)
        longest = 0
        for i in range(n):
            while j < n and s[j] not in unique_chars:
                unique_chars.add(s[j])
                j += 1
            longest = max(longest, j - i)
            unique_chars.remove(s[j])
        return longest


class Solution2:
    def lengthOfLongestSubstring(self, s):
        st = {}
        i, ans = 0, 0
        for j in range(len(s)):
            if s[j] in st:
                i = max(st[s[j]], i)
            ans = max(ans, j - i + 1)
            st[s[j]] = j + 1
        return ans


#Solution3:滑动窗口法
class Solution3:
    def lengthOfLongestSubstring(self, s):
        table = dict()
        slow = 0
        fast = 0
        ret_max = 0

        while slow < len(s) and fast < len(s):
            if s[fast] in table:
                del(table[s[slow]])
                slow += 1
            else:
                table[s[fast]] = fast
                fast += 1
                ret_max = max(ret_max, fast - slow)
        return ret_max