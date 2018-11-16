'''
2.两数相加
给定两个非空链表来表示两个非负整数。位数按照逆序方式存储，它们的每个节点只存储单个数字。
将两数相加返回一个新的链表。
你可以假设除了数字 0 之外，这两个数字都不会以零开头。

示例：

输入：(2 -> 4 -> 3) + (5 -> 6 -> 4)
输出：7 -> 0 -> 8
原因：342 + 465 = 807

'''
#这道题目涉及链表的相关知识点，未搞明白
#在有空闲时间继续多看几遍，争取搞懂

class ListNode:
    def __init__(self, x):
        self.val = x
        self.next = None


class Solution:
    def addTwoNmbers(self, l1, l2):
        head = ListNode(0)
        ptr = head
        carry = 0
        while True:
            if l1 != None:
                carry += l1.val
                l1 = l1.next
            if l2 != None:
                carry += l2.val
                l2 = le.next
            ptr.val = carry % 10
            carry /=10

            if l1 != None or l2 != None or carry != 0:
                ptr.next = ListNode(0)
                ptr = ptr.next
            else:
                break
        return head