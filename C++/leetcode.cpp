#include <iostream>
#include <vector>
#include <stack>
#include <string>
#include <math.h>
#include <unordered_map>
#include <sstream>

using namespace std;

struct ListNode {
    int val;
    ListNode *next;

    ListNode() : val(0), next(nullptr) {}

    ListNode(int x) : val(x), next(nullptr) {}

    ListNode(int x, ListNode *next) : val(x), next(next) {}
};

struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;

    TreeNode() : val(0), left(nullptr), right(nullptr) {}

    TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}

    TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}
};


class Solution {
public:
    static vector<int> twoSum(vector<int> &nums, int target) {
        int sum = 0;
        vector<int> ans;
        for (int i = 0; i < nums.size(); i++) {
            for (int j = i + 1; j < nums.size(); j++) {
                if (nums[i] + nums[j] == target) {
                    ans.push_back(i);
                    ans.push_back(j);
                    break;
                }
            }
        }
        return ans;
    }

    static bool isPalindrome(const int x) {
        string ss = to_string(x);
        string s;
        std::reverse(ss.begin(), ss.end());
        cout << ss << '\n';
        s = ss;
        std::reverse(ss.begin(), ss.end());
        if (x < 0) {
            s += '-';
        }
        if (s == ss) {
            return true;
        } else {
            return false;
        }
    }

    string longestCommonPrefix(vector<string> &strs) {
        string ans = "";
        sort(strs.begin(), strs.end());
        string v = strs[0];
        string v2 = strs[strs.size() - 1];
        for (int i = 0; i < min(v.size(), v2.size()); i++) {
            if (v[i] != v2[i]) {
                return ans;
            }
            ans += v[i];
        }
        return ans;
    }

    static bool isValid(string s) {
        stack<char> st;
        for (int i = 0; i < s.size(); i++) {
            if (s[i] == '(' || s[i] == '[' || s[i] == '{') {
                st.push(s[i]);
            } else if (st.empty() || (s[i] != ')' && st.top() == '(') || (s[i] != ']' && st.top() == '[') ||
                       (s[i] != '}' && st.top() == '{')) {
                return false;
            } else {
                st.pop();
            }
        }
        return st.empty();
    }

    ListNode *mergeTwoLists(ListNode *list1, ListNode *list2) {
        if (list1 == NULL) {
            return list2;
        }
        if (list2 == NULL) {
            return list1;
        }
        if (list1->val <= list2->val) {
            list1->next = mergeTwoLists(list1->next, list2);
            return list1;
        } else {
            list2->next = mergeTwoLists(list1, list2->next);
            return list2;
        }
    }

    int removeDuplicates(vector<int> &nums) {
        int i = 0;
        for (int j = 1; j < nums.size(); j++) {
            if (nums[i] != nums[j]) {
                nums[i + 1] = nums[j];
                i++;
            }
        }
        return i + 1;
    }

    int removeElement(vector<int> &nums, int val) {
        int ans = 0;
        for (int i = 0; i < nums.size(); i++) {
            if (nums[i] == val) {
                nums[ans++] = nums[i];
            }
        }
        return ans;
    }

    int strStr(string haystack, string needle) {
        int a = haystack.find(needle);
        if (a != string::npos) {
            return a;
        } else {
            return -1;
        }
    }

    int searchInsert(vector<int> &nums, int target) {
        int l = 0;
        int r = nums.size();
        int m;
        if (target > nums[r - 1]) {
            return r;
        }
        while (l <= r) {
            m = (l + r) / 2;
            if (nums[m] == target) {
                return m;
            }
            if (target < nums[m]) {
                r = m - 1;
            } else {
                l = m + 1;
            }
        }
        return l;
    }

    int lengthOfLastWord(string s) {
        int cnt = s.size() - 1;
        while (s[cnt] == ' ') {
            cnt--;
        }
        int c = 0;
        for (; cnt >= 0; cnt--) {
            if (s[cnt] == ' ') {
                return c;
            }
            c++;
        }
        return c;
    }
//    vector<int> plusOne(vector<int>& digits) {
//
//    }

    static int mySqrt(int x) {
        return floor(exp(0.5 * log(x)));
    }

    static string addBinary(string a, string b) {
        string ans;
        int carry = 0;
        int i = a.length() - 1;
        int j = b.length() - 1;

        while (i >= 0 || j >= 0 || carry) {
            if (i >= 0)
                carry += a[i--] - '0';
            if (j >= 0)
                carry += b[j--] - '0';
            ans += carry % 2 + '0';
            carry /= 2;
        }

        std::reverse(ans.begin(), ans.end());
        return ans;
    }

    int singleNumber(vector<int> &nums) {
        sort(nums.begin(), nums.end());
        for (int i = 1; i < nums.size(); i += 2) {
            if (nums[i] != nums[i - 1]) {
                return nums[i - 1];
            }
        }
        return nums[nums.size() - 1];
    }

    vector<string> summaryRanges(vector<int> &nums) {
        vector<string> result;
        int n = nums.size();
        if (n == 0) return result;
        int a = nums[0];
        for (int i = 0; i < n; ++i) {
            if (i == n - 1 || nums[i] + 1 != nums[i + 1]) {
                if (nums[i] != a) {
                    result.push_back(to_string(a) + "->" + to_string(nums[i]));
                } else {
                    result.push_back(to_string(a));
                }
                if (i != n - 1) {
                    a = nums[i + 1];
                }
            }
        }
        return result;
    }

    ListNode *reverseList(ListNode *head) {
        ListNode *prev = NULL;
        ListNode *curr = head;
        while (curr != NULL) {
            ListNode *forward = curr->next;
            curr->next = prev;
            prev = curr;
            curr = forward;
        }
        return prev;
    }

    void moveZeroes(vector<int> &nums) {
        int cnt = 0;
        int n = nums.size();
        for (int i = 0; i < n; i++) {
            if (nums[i] != 0) {
                nums[cnt++] = nums[i];
            }
        }
        while (cnt < n) {
            nums[cnt++] = 0;
        }

    }

    bool isPalindrome(string s) {
        string str;
        for (int i = 0; i < s.size(); i++) {
            if (isupper(s[i])) {
                s[i] = tolower(s[i]);
            }
            if (isalnum(s[i])) {
                str.push_back(s[i]);
            }
        }
        if (string(str.rbegin(), str.rend()) == str) {
            return true;
        } else {
            return false;
        }
    }

    static int missingNumber(vector<int> &nums) {
        int n = nums.size() + 1;
        int t = (n * (n - 1)) / 2;
        for (auto x: nums) {
            t -= x;
        }
        return t;
    }

    bool check(TreeNode *left, TreeNode *right) {
        if (left == NULL || right == NULL) {
            return left == right;
        }
        if (left->val != right->val) {
            return false;
        }
        return check(left->left, right->right) && check(left->right, right->left);

    }

    bool isSymmetric(TreeNode *root) {
        if (root == NULL) {
            return true;
        }
        return check(root->left, root->right);
    }


    vector<int> intersect(vector<int> &nums1, vector<int> &nums2) {
        sort(nums1.begin(), nums1.end());
        sort(nums2.begin(), nums2.end());
        int n = nums1.size(), m = nums2.size();
        vector<int> ans;
        int i = 0, j = 0;
        while (i < n && j < m) {
            if (nums1[i] == nums2[j]) {
                ans.push_back(nums1[i]);
                i++, j++;
            } else if (nums1[i] < nums2[j]) {
                i++;
            } else {
                j++;
            }
        }
        return ans;
    }

    vector<int> sortedSquares(vector<int> &nums) {
        for (int i = 0; i < nums.size(); i++) {
            nums[i] *= nums[i];
        }
        sort(nums.begin(), nums.end());
        return nums;
    }

    bool isSubsequence(string s, string t) {
        int n = s.length();
        int m = t.length();
        int i = 0, j = 0;
        while (i < n && j < m) {
            if (s[i] == t[j]) {
                i++;
            }
            j++;
        }
        return i == n;
    }

    void solve(int n, int open, int close, string s, vector<string> ans) {
        if (open == close && close == n) {
            ans.push_back(s);
            return;
        }
        if (open < n) {
            solve(n, open + 1, close, s + "(", ans);
        }
        if (close < n) {
            solve(n, open, close + 1, s + ")", ans);
        }

    }

    vector<string> generateParenthesis(int n) {
        vector<string> ans;
        solve(n, 0, 0, "", ans);
        return ans;
    }

    int longestSubarray(vector<int> &nums) {
        int res = 0, cnt = 0, i = 0, j = 0;
        while (j < nums.size()) {
            if (nums[j] == 0) {
                cnt++;
            }
            while (cnt > 1) {
                if (nums[i] == 0) {
                    cnt--;
                }
                i++;
            }
            res = max(res, j - i);
            j++;
        }
        return res;
    }

    int longestOnes(vector<int> &nums, int k) {
        int start = 0;
        int i = 0;
        while (start < nums.size()) {
            if (nums[start] == 0) {
                k--;
            }
            if (k < 0) {
                if (nums[i] == 0) {
                    k++;
                }
                i++;
            }
            start++;
        }
        return start - i;
    }

    ListNode *addTwoNumbers(ListNode *l1, ListNode *l2) {
        ListNode *result = new ListNode(0);
        ListNode *cn1 = l1;
        ListNode *cn2 = l2;
        ListNode *curr = result;
        int if_zero = 0;
        while (cn1 != NULL || cn2 != NULL) {
            int val1 = (cn1 == NULL) ? 0 : cn1->val;
            int val2 = (cn2 == NULL) ? 0 : cn2->val;
            int sum = val1 + val2 + if_zero;
            if_zero = sum / 10;
            curr->next = new ListNode(sum % 10);
            curr = curr->next;
            if (cn2 != NULL) {
                cn2 = cn2->next;
            }
            if (cn1 != NULL) {
                cn1 = cn1->next;
            }
        }
        if (if_zero > 0) {
            curr->next = new ListNode(if_zero);
        }
        return result->next;
    }

    int begin;

    void check(int l, int r, string &s, int &max, string &ans) {
        while (l >= 0 && r < s.size() && s[l] == s[r]) {
            l--;
            r++;
        }
        l++;
        r--;

        int a = r - l + 1;
        if (a > max) {
            ans = s.substr(l, a);
            max = a;
        }
    }

    string longestPalindrome(string s) {
        int n = s.length();
        int max = 0;
        string ans = "";
        for (int i = 0; i < n; i++) {
            check(i, i, s, max, ans);
            if (i == n - 1) {
                break;
            }
            check(i, i + 1, s, max, ans);
            if (max == n) {
                return s;
            }
        }
        return ans;
    }


    void rotate(vector<vector<int>> &matrix) {
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j <= i; j++) {
                swap(matrix[i][j], matrix[j][i]);
            }

        }
        for (int i = 0; i < matrix.size(); i++) {
            std::reverse(matrix[i].begin(), matrix[i].end());
        }
    }


    vector<vector<string>> groupAnagrams(vector<string> &strs) {
        unordered_map<string, vector<string>> ans;
        for (auto str: strs) {
            string s = str;
            sort(s.begin(), s.end());
            ans[s].push_back(str);
        }
        vector<vector<string>> res;
        for (auto str: ans) {
            res.push_back(str.second);
        }
        return res;
    }

    bool chck(TreeNode *root, int l, int r) {
        if (root->val < r && root->val > l) {
            return chck(root->left, l, root->val) && chck(root->right, root->val, r);
        } else {
            return false;
        }
    }

    bool isValidBST(TreeNode *root) {
        int min = -1 * INT_MAX, max = INT_MAX;
        return chck(root, min, max);
    }

    ListNode *deleteDuplicates(ListNode *head) {
        ListNode *ans = head;
        while (ans->next && ans) {
            if (ans->val == ans->next->val) {
                ans->next = ans->next->next;
                continue;
            }
            ans = ans->next;
        }
        return head;
    }

    int majorityElement(vector<int> &nums) {
        int n = nums.size();
        unordered_map<int, int> m;
        for (int i = 0; i < n; i++) {
            m[nums[i]]++;
        }
        n /= 2;
        for (auto x: m) {
            if (x.second > n) {
                return x.first;
            }
        }
        return 0;
    }

    vector<vector<int>> generate(int numRows) {
        vector<vector<int>> ans;
        for (int i = 0; i < numRows; i++) {
            vector<int> row(i + 1, i);
            for (int j = 1; j < i; j++) {
                row[j] = ans[i - 1][j] + ans[i - 1][j - 1];
            }
            ans.push_back(row);
        }
        return ans;
    }

    vector<int> getRow(int rowIndex) {
        vector<int> row(rowIndex + 1, 1);
        for (int i = 0; i <= rowIndex; i++) {
            for (int j = i - 1; j > 0; j--) {
                row[j] = row[j] + row[j - 1];
            }
        }
        return row;
    }

    int maxProfit(vector<int> &prices) {
        int n = prices.size();
        int buy = prices[0], sell = 0;
        for (int i = 0; i < n; i++) {
            buy = min(buy, prices[i]);
            sell = max(sell, prices[i] - buy);
        }
        return sell;
    }

    int reverse(int x) {
        long long xx = static_cast<long long>(x);
        bool flag = false;
        if (xx < 0) {
            flag = true;
            xx *= -1;
        }
        string s = to_string(xx);
        if (s[s.size()] == '0') {
            s.pop_back();
        }
        std::reverse(s.begin(), s.end());
        if (flag) {
            s = '-' + s;
        }
        if (stoll(s) > INT_MAX || stoll(s) < INT_MIN) {
            return 0;
        }
        return stoll(s);
    }

    int myAtoi(string s) {
        stringstream st(s);
        int i = 0;
        st >> i;
        return i;
    }

    int lengthOfLongestSubstring(string s) {
        int n = s.length();
        int cnt[128] = {};
        int ans = 0;
        for (int j = 0, i = 0; j < n; j++) {
            cnt[s[j]]++;
            while (cnt[s[j]] > 1) {
                cnt[s[i++]]--;
            }
            ans = max(ans, j - i + 1);
        }
        return ans;
    }

    string sortSentence(string s) {
        string str;
        vector<pair<int, string>> vec;
        istringstream ss(s);
        string word;
        while (ss >> word) {
            int ind = word[word.size() - 1];
            word.pop_back();
            vec.push_back({ind, word});
        }
        sort(vec.begin(), vec.end());
        for (auto it: vec) {
            str += it.second;
            str += ' ';
        }
        str.pop_back();
        return str;
    }

    int findMin(vector<int> &nums) {
        int n = nums.size();
        int l = 0, r = n - 1;
        while (l < r) {
            int m = (r - l) / 2;
            if (nums[m] >= nums[r]) {
                l = m + 1;
            } else {
                r = m;
            }
        }
        return nums[l];
    }

    string addSpaces(string s, vector<int> &spaces) {
        int n = spaces.size(), m = s.size(), i = 0, j = 0;
        string ans = "";
        while (i < m) {
            if (j < n && i == spaces[j]) {
                ans += " ";
                j++;
            }
            ans += s[i];
            i++;
        }
        return ans;
    }
//    vector<int> twoSum(vector<int>& numbers, int target) {
//        int n =numbers.size();
//        for(int i = 0, j = n-1; i < j; ){
//            int sum = numbers[i] + numbers[j];
//            if(sum == target) return {i+1, j +1};
//            else if(sum < target) i++;
//            else j--;
//        }
//        return {};
//    }

    bool isPossibleDivide(vector<int> &nums, int k) {
        if (nums.size() % k != 0) {
            return false;
        }
        unordered_map<int, int> mapik;
        for (int i = 0; i < nums.size(); i++) {
            mapik[nums[i]]++;
        }
        sort(nums.begin(), nums.end());
        for (auto num: nums) {
            if (mapik[num] > 0) {
                for (int i = num + 1; i < num + k; i++) {
                    if (mapik[i] == 0) {
                        return false;
                    }
                    mapik[i]--;
                }
                mapik[num]--;

            }
        }
        return true;
    }

    vector<int> findKDistantIndices(vector<int> &nums, int key, int k) {
        vector<int> ans;
        for (int i = 0; i < nums.size(); i++) {
            int a = -1;
            for (int j = i; j < nums.size(); j++) {
                if (abs(i - j) <= k && nums[j] == k) {
                    a = i;
                }

            }
            if (a != -1) {
                ans.push_back(a);
            }
        }
        return ans;
    }

    static int numOfPairs(vector<string> &nums, string target) {
        int ans = 0;
        long long k = stoll(target);
        for (int i = 0; i < nums.size(); i++) {
            for (int j = 0; j < nums.size(); j++) {
                if (i != j && stoll(nums[i] + nums[j]) == k) {
                    ans++;
                }
            }
        }
        return ans;
    }

    vector<int> findDisappearedNumbers(vector<int> &nums) {
        unordered_map<int, int> mapchik;
        vector<int> ans;
        for (auto x: nums) {
            mapchik[x]++;
        }
        for (int i = 1; i <= nums.size(); i++) {
            if (mapchik[i] == 0) {
                ans.push_back(i);
            }
        }
        return ans;
    }


    int calculate(string s) {
        istringstream in('+' + s + '+');
        int ans = 0, cnt = 0, n;
        char znak;
        while (in >> znak) {
            if (znak == '-' || znak == '+') {
                ans += cnt;
                in >> cnt;
                cnt *= 44 - znak;
            } else {
                in >> n;
                if (znak == '*') {
                    cnt *= n;
                } else {
                    cnt /= n;
                }
            }
        }
        return ans;
    }

    int titleToNumber(string columnTitle) {
        int ans = 0;
        for (auto s: columnTitle) {
            int d = s - 'A' + 1;
            ans = ans * 26 + d;
        }
        return ans;
    }

    int gcd_help(int a, int b) {
        if (a == 0) {
            return b;
        }
        return gcd_help(b % a, a);
    }

    int nthMagicalNumber(int n, int a, int b) {
        long long int l = min(a, b);
        long long r = n * l;
        long long lcm = (a * b) / gcd_help(a, b);
        while (l < r) {
            long long m = l + (r - l) / 2;
            long long fac = m / a + m / b + m / lcm;
            if (fac < n) {
                l = m + 1;
            } else {
                r = m;
            }
        }
        int mod = 1e9 + 7;
        return l % mod;
    }

    int minIncrementForUnique(vector<int> &nums) {
        int ans = 0;
        sort(nums.begin(), nums.end());
        for (int i = 1; i < nums.size(); i++) {
            if (nums[i] <= nums[i - 1]) {
                ans += abs(nums[i] - nums[i - 1]) + 1;
                nums[i] += abs(nums[i] - nums[i - 1]) + 1;
            }
        }
        return ans;
    }

    string orderlyQueue(string s, int k) {
        string ans = s;
        if (k > 1) {
            sort(s.begin(), s.end());
            return s;
        } else {
            string temp = "";
            temp = s + s;
            for (int i = 0; i + s.length() <= temp.length(); i++) {
                ans = min(ans, temp.substr(i, s.length()));
            }
        }
        return ans;
    }

    double minAreaFreeRect(vector<vector<int>> &points) {
        if (points.size() < 4) {
            return 0;
        }
        double ans = 0.0;
        int x0, y0;
        int x1, y1;
        int x2, y2;
        int x3, y3;
        int Lx1, Lx2;
        int Ly1, Ly2;
        for (int i = 0; i < points.size() - 3; i++) {
            x0 = points[i][0];
            y0 = points[i][1];
            for (int j = i + 1; j < points.size(); j++) {
                x1 = points[j][0];
                y1 = points[j][1];
                for (int k = j + 1; k < points.size(); k++) {
                    x2 = points[k][0];
                    y2 = points[k][1];
                    Lx1 = x1 - x0;
                    Ly1 = y1 - y0;

                    Lx2 = x2 - x0;
                    Ly2 = y2 - y0;
                    int dotProd = Lx1 * Lx2 + Ly1 * Ly2;
                    if (dotProd != 0) {
                        continue;
                    }
                    bool fl = true;
                    for (int n = 0; n < points.size(); n++) {
                        x3 = points[n][0];
                        y3 = points[n][1];
                        if ((x3 == x0 + Lx1 + Lx2) && (y3 == y0 + Ly1 + Ly2)) {
                            fl = false;
                            break;
                        }
                    }
                    if (fl) {
                        continue;
                    }
                    double area = (double) abs(Lx1 * Ly2 - Ly1 * Lx2);
                    if (area == 0.0) {
                        area = ans;
                    } else {
                        area = area < ans ? area : ans;
                    }


                }
            }
        }
        return ans;
    }

    int minCostToMoveChips(vector<int> &position) {
        int odd = 0;
        int ev = 0;
        for (int i = 0; i < position.size(); i++) {
            if (position[i] % 2 == 0) {
                ev++;
            } else {
                odd++;
            }

        }
        return min(ev, odd);
    }

    int maxArea(vector<int> &height) {
        int ans = 0;
        int l = 0, r = height.size() - 1;
        while (l <= r) {
            int minn = 10329084234;
            minn = min(min(height[l],height[r]), minn);
            ans = max(minn*(r-l), ans);
            if(height[l] > height[r]){
                r--;
            } else {
                l++;
            }
        }
        return ans;
    }


};


int main() {
    cout << '\n';
    bool flag = Solution::isPalindrome(10);
    bool fl = Solution::isValid("{}[]()");
    cout << fl << '\n';
    cout << flag;
    cout << Solution::mySqrt(4) << '\n';
    cout << Solution::addBinary("11", "1") << '\n';
    vector<int> c = {0};
    vector<string> nn = {"777", "7", "77", "77"};
    cout << Solution::numOfPairs(nn, "7777") << '\n';
    cout << Solution::missingNumber(c) << '\n';
}
