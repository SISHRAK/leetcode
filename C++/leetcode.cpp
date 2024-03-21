#include <iostream>
#include <vector>
#include <stack>
#include <string>
#include <map>
#include <math.h>
#include <unordered_map>
#include <queue>
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

    [[maybe_unused]] void moveZeroes(vector<int> &nums) {
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
            minn = min(min(height[l], height[r]), minn);
            ans = max(minn * (r - l), ans);
            if (height[l] > height[r]) {
                r--;
            } else {
                l++;
            }
        }
        return ans;
    }

    int longestSubarray_1(vector<int> &nums) {
        int ans = 0;
        int maxx_ = 0;
        int i = 0, j = 0;
        int ct0 = 0;
        while (i < nums.size()) {
            if (nums[i] == 0) {
                ct0++;
                while (ct0 > 1) {
                    if (nums[j] == 0) {
                        ct0--;
                    } else {
                        ans--;
                    }
                    j++;
                }
            } else {
                ans++;
                maxx_ = max(ans, maxx_);
            }
            i++;
        }
        if (maxx_ == nums.size()) {
            return --maxx_;
        } else {
            return maxx_;
        }
    }

    vector<string> summaryRangess(vector<int> &nums) {
        vector<string> ans;
        int n = nums.size();
        if (n == 0) {
            return ans;
        }
        int a = nums[0];
        for (int i = 1; i < n; i++) {
            if (i == n - 1 || nums[i] + 1 != nums[i + 1]) {
                if (nums[i] != a) {
                    ans.push_back(to_string(a) + "->" + to_string(nums[i]));
                } else {
                    ans.push_back(to_string(a));
                }
            }
            if (i != n - 1) {
                a = nums[i + 1];
            }
        }
        return ans;
    }

    static int compress(vector<char> &chars) {
        int ans = 0;
        for (int i = 0; i < chars.size(); i++) {
            int cnt = 0;
            char let = chars[i];
            while (i < chars.size() && chars[i] == let) {
                ++cnt;
                ++i;
            }
            chars[ans++] = let;
            if (cnt > 1) {
                for (auto key: to_string(cnt)) {
                    chars[ans++] = key;
                }
            }

        }
        return chars.size();
    }

    bool isPalindromee(string s) {
        string s1 = "";
        for (char &i: s) {
            if (isupper(i)) {
                i = tolower(i);
            }
            if (isalnum(i)) {
                s1.push_back(i);
            }
        }
        if (string(s1.rbegin(), s1.rend()) == s1) {
            return true;
        } else {
            return false;
        }
    }

    int subarraySum(vector<int> &nums, int k) {
        map<int, int> anss;
        int ans = 0;
        int sum = 0;
        anss[sum] = 1;
        for (auto e: nums) {
            sum += e;
            ans += anss[sum - k];
            anss[sum]++;
        }
        return ans;
    }

    void moveZeroess(vector<int> &nums) {
        int cnt = 0;
        for (int i = 0; i < nums.size(); i++) {
            if (nums[i] != 0) {
                nums[cnt++] = nums[i];
            }
        }
        while (cnt < nums.size()) {
            nums[cnt++] = 0;
        }
    }

    vector<vector<string>> groupAnagramss(vector<string> &strs) {
        map<string, vector<string>> mapp;
        for (auto s: strs) {
            string str = s;
            sort(str.begin(), str.end());
            mapp[str].push_back(s);
        }
        vector<vector<string>> ans;
        for (auto k: mapp) {
            ans.push_back(k.second);
        }
        return ans;
    }

    void solve2(int n, int open, int closed, string s, vector<string> &a) {
        if (open == n && closed == n) {
            a.push_back(s);
            return;
        }
        if (open < n) {
            solve2(n, open + 1, closed, s + "(", a);
        }
        if (closed < open) {
            solve2(n, open, closed + 1, s + ")", a);
        }

    }

    vector<string> generateParenthesiss(int n) {
        vector<string> a;
        solve2(n, 0, 0, "", a);
        return a;
    }

    ListNode *reverseLists(ListNode *head) {
        ListNode *prev = NULL;
        ListNode *curr = head;
        while (curr != NULL) {
            ListNode *forw = curr->next;
            curr->next = prev;
            prev = curr;
            curr = forw;
        }
        return prev;
    }

    bool checkInclusion(string s1, string s2) {
        int m1 = s1.length();
        int i = 0;
        string str = s2.substr(i, m1);
        sort(s1.begin(), s2.end());
        while (str.length() == m1) {
            str = s2.substr(i, m1);
            sort(str.begin(), str.end());
            if (str == s1) {
                return true;
            }
            i++;
        }
        return false;
    }

    bool isValidd(string s) {
        stack<char> st;
        for (int i = 0; i < s.size(); i++) {
            if (s[i] == '{' || s[i] == '[' || s[i] == '(') {
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

//    ListNode *reverseList(ListNode *head) {
//        ListNode *prev = NULL;
//        ListNode *curr = head;
//        while (curr != NULL) {
//            ListNode *forward = curr->next;
//            curr->next = prev;
//            prev = curr;
//            curr = forward;
//        }
//        return prev;
//    }
//1 1 1
};

class RecentCounter {
public:
    queue<int> a;

    RecentCounter() {
    }

    int ping(int t) {

        a.push(t);
        while (a.front() < t - 3000) {
            a.pop();
        }
        return a.size();
    }

    int maxDistToClosest(vector<int> &seats) {
        int ans = 0;
        int l, r;
        for (int i = 0; i < seats.size(); i++) {
            l = i;
            r = i;
            while (l >= 0 && seats[l] != 1) l--;
            while (r < seats.size() && seats[r] != 1) r++;
            if (l < 0) {
                l = r;
            }
            if (r == seats.size()) {
                r = l;
            }
            ans = max(ans, min(abs(l - i), abs(r - i)));

        }
        return ans;
    }

    vector<vector<int>> merge(vector<vector<int>> &intervals) {
        int n = intervals.size();
        if (n <= 1) {
            return intervals;
        }
        sort(intervals.begin(), intervals.end());
        vector<vector<int>> ans;
        ans.push_back(intervals[0]);
        for (int i = 1; i < n; i++) {
            if (intervals[i][0] <= ans.back()[1]) {
                ans.back()[1] = max(ans.back()[1], intervals[i][1]);
            } else {
                ans.push_back(intervals[i]);
            }
        }
        return ans;
    }

    int trap(vector<int> &height) {
        int ans = 0;
        int n = height.size();
        vector<int> pref(n);
        vector<int> sub(n);
        pref[0] = height[0];
        for (int i = 1; i < n; i++) {
            pref[i] = max(height[i], pref[i - 1]);
        }
        sub[n - 1] = height[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            sub[i] = max(height[i], sub[i + 1]);
        }
        for (int i = 0; i < n; i++) {
            ans += min(pref[i], sub[i]) - height[i];
        }
        return ans;
    }

    vector<int> twoSum(vector<int> &nums, int target) {
        unordered_map<int, int> aa;
        int num;
        for (int i = 0; i < nums.size(); i++) {
            num = nums[i];
            if (aa.find(target - num) != aa.end()) {
                return {i, aa[target - num]};
            }
            aa[num] = i;
        }
        return {};
    }

    vector<int> findAnagrams(string s, string p) {
        vector<int> ans;
        if (s.size() < p.size()) {
            return ans;
        }
        vector<int> a(26);
        vector<int> b(26);
        for (auto c: p) {
            b[c - 'a']++;
        }
        int i = 0, j = 0;
        while (j < s.size()) {
            a[s[j] - 'a']++;
            if (j - i + 1 == p.size()) {
                if (a == b)
                    ans.push_back(i);
            }
            if (p.size() > j - i + 1) {
                j++;
            } else {
                a[s[i] - 'a']--;
                i++;
                j++;
            }
        }
        return ans;
    }

    int rand7();

    int rand10() {
        int n = rand7(), m = 7;
        while (n > 5) {
            m = n - 5;
            n = rand7();
        }
        while (m == 7) m = rand7();
        return (m % 2 ? 5 : 0) + n;
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

    bool chck(TreeNode *root, long long l, long long r) {
        if (root == NULL) {
            return true;
        }
        if (root->val < r && root->val > l) {
            return chck(root->left, l, root->val) && chck(root->right, root->val, r);
        } else {
            return false;
        }
    }

    bool isValidBST(TreeNode *root) {
        long long int min = -100000000000, max = 10000000000;
        return chck(root, min, max);
    }

    void dfs(int i, int j, vector<vector<char>> &grid) {
        grid[i][j] = 0;
        int n = grid.size();
        int m = grid[0].size();

        if (i - 1 >= 0 && grid[i - 1][j] == '1') dfs(i - 1, j, grid);
        if (i + 1 < n && grid[i + 1][j] == '1') dfs(i + 1, j, grid);
        if (j - 1 >= 0 && grid[i][j - 1] == '1') dfs(i, j - 1, grid);
        if (j + 1 < m && grid[i][j + 1] == '1') dfs(i, j + 1, grid);

    }

    int numIslands(vector<vector<char>> &grid) {
        int n = grid.size();
        int m = grid[0].size();
        int ans = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (grid[i][j] == '1') {
                    ans++;
                    dfs(i, j, grid);
                }
            }
        }
        return ans;
    }

    int maxPower(string s) {
        int cnt = 1, ans = 1;
        char buf = s[0];
        for (int i = 1; i < s.size(); i++) {
            char cur = s[i];
            if (cur == buf) cnt++;
            else {
                buf = cur;
                cnt = 1;
            }
            ans = max(ans, cnt);
        }
        return ans;
    }

    vector<vector<int>> intervalIntersection(vector<vector<int>> &firstList, vector<vector<int>> &secondList) {
        vector<vector<int>> ans;
        for (int i = 0, j = 0; i < firstList.size() && j < secondList.size();) {
            if (firstList[i][1] < secondList[j][0]) {
                i++;
            } else if (secondList[j][1] < firstList[i][0]) {
                j++;
            } else {
                ans.push_back({max(firstList[i][0], secondList[j][0]), min(firstList[i][1], secondList[j][1])});
                if (firstList[i][1] < secondList[j][1]) {
                    i++;
                } else {
                    j++;
                }
            }
        }
        return ans;
    }

    ListNode *addTwoNumbers(ListNode *l1, ListNode *l2) {
        auto *dum = new ListNode();
        ListNode *temp = dum;
        int curry = 0;
        while (l1 != NULL || l2 != NULL || curry) {
            int sum = 0;
            if (l1 != NULL) {
                sum += l1->val;
                l1 = l1->next;
            }
            if (l2 != NULL) {
                sum += l2->val;
                l2 = l2->next;
            }
            sum += curry;
            curry = sum / 10;
            ListNode *newnode = new ListNode(sum % 10);
            temp->next = newnode;
            temp = temp->next;
        }
        return dum->next;
    }

    void merge(vector<int> &nums1, int m, vector<int> &nums2, int n) {
        for (int j = 0, i = m; j < n; j++) {
            nums1[i] = nums2[j];
            i++;
        }
        sort(nums1.begin(), nums1.end());
    }

    ListNode *mergeTwoLists(ListNode *list1, ListNode *list2) {
        if (list1 == NULL) {
            return list2;
        }
        if (list2 == NULL) {
            return list1;
        }
        if (list1->val < list2->val) {
            list1->next = mergeTwoLists(list1->next, list2);
            return list1;
        } else {
            list2->next = mergeTwoLists(list1, list2->next);
            return list2;
        }
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

    string longestPalindrome(string s) {
        if (s.empty()) {
            return "";
        }
        int start = 0;
        int end = 0;
        for (int i = 0; i < s.size(); i++) {
            int odd = ExpandAroundCenter(s, i, i);
            int even = ExpandAroundCenter(s, i, i + 1);
            int max_len = max(odd, even);
            if (max_len > end - start) {
                start = i - (max_len - 1) / 2;
                end = i + max_len / 2;
            }
        }
        return s.substr(start, end - start + 1);
    }

    int ExpandAroundCenter(string s, int left, int right) {
        while (left >= 0 && right < s.size() && s[left] == s[right]) {
            left--;
            right++;
        }
        return right - left - 1;
    }

    int numJewelsInStones(string jewels, string stones) {
        int count = 0;
        for (char jewel: jewels) {
            for (char stone: stones) {
                if (jewel == stone) {
                    count++;
                }
            }
        }
        return count;
    }

    TreeNode *lowestCommonAncestor(TreeNode *root, TreeNode *p, TreeNode *q) {
        if (root == NULL) {
            return NULL;
        }
        if (root == p || root == q) {
            return root;
        }
        TreeNode *l = lowestCommonAncestor(root->left, p, q);
        TreeNode *r = lowestCommonAncestor(root->right, p, q);
        if (l && r) {
            return root;
        }
        if (l) {
            return l;
        }
        if (r) {
            return r;
        }
        return NULL;
    }

    vector<int> intersect(vector<int> &nums1, vector<int> &nums2) {
        sort(nums1.begin(), nums1.end());
        sort(nums2.begin(), nums2.end());
        int i = 0, j = 0;
        vector<int> ans;
        while (i < nums1.size() && j < nums2.size()) {
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

    int missingNumber(vector<int> &nums) {
        int n = nums.size() + 1;
        int t = (n * (n - 1)) / 2;
        for (auto x: nums) {
            t -= x;
        }
        return t;
    }

    int evalRPN(vector<string> &tokens) {
        stack<string> st;
        for (auto s: tokens) {
            if (isOperator(s)) {
                int res = 0;
                int el2 = stoi(st.top());
                st.pop();
                int el1 = stoi(st.top());
                st.pop();
                if (s == "+") {
                    res = el1 + el2;
                }
                if (s == "*") {
                    res = el1 * el2;
                }
                if (s == "-") {
                    res = el1 - el2;
                }
                if (s == "/") {
                    res = el1 / el2;
                }
                st.push(to_string(res));
            } else {
                st.push(s);
            }
        }
        return stoi(st.top());
    }

    bool isOperator(string s) {
        return (s == "-" || s == "+" || s == "/" || s == "*");
    }

    double findMedianSortedArrays(vector<int> &nums1, vector<int> &nums2) {
        int n = nums1.size();
        int m = nums2.size();
        int i = 0, j = 0, m1 = 0, m2 = 0;
        for (int count = 0; count <= (n + m) / 2; count++) {
            m2 = m1;
            if (i != n && j != m) {
                if (nums1[i] > nums2[j]) {
                    m1 = nums2[j++];
                } else {
                    m1 = nums1[i++];
                }
            } else if (i < n) {
                m1 = nums1[i++];
            } else {
                m1 = nums2[j++];
            }
        }
        if ((n + m) % 2 == 1) {
            return static_cast<double>(m1);
        } else {
            double ans = static_cast<double>(m1) + static_cast<double>(m2);
            return ans / 2.0;
        }
    }

    string simplifyPath(string path) {
        stack<string> st;
        string ans = "";
        for (int i = 0; i < path.size(); i++) {
            if (path[i] == '/') {
                continue;
            }
            string temp;
            while (i < path.size() && path[i] != '/') {
                temp += path[i];
                i++;
            }
            if (temp == ".") {
                continue;
            } else if (temp == "..") {
                if (!st.empty()) {
                    st.pop();
                }
            } else {
                st.push(temp);
            }
        }
        while (!st.empty()) {
            ans = '/' + st.top() + ans;
            st.pop();
        }
        if (ans.size() == 0) {
            return "/";
        }
        return ans;
    }

    bool isSubsequence(string s, string t) {
        int pos = 0;
        for (char i: t) {
            if (i == s[pos]) {
                pos++;
            }
        }
        return pos == s.size();
    }

    vector<int> sortedSquares(vector<int> &nums) {
        int n = nums.size();
        vector<int> ans(n);
        for (int i = 0; i < n; i++) {
            nums[i] = nums[i] * nums[i];
        }
        int left = 0;
        int right = n - 1;
        for (int i = n - 1; i >= 0; i--) {
            if (nums[left] > nums[right]) {
                ans[i] = nums[left];
                left++;
            } else {
                ans[i] = nums[right];
                right--;
            }
        }
        return ans;
    }

    ListNode *removeNthFromEnd(ListNode *head, int n) {
        ListNode *dum = new ListNode(0);
        dum->next = head;
        ListNode *first = dum;
        ListNode *second = dum;

        for (int i = 0; i <= n; i++) {
            first = first->next;
        }
        while (first != NULL) {
            first = first->next;
            second = second->next;
        }
        ListNode *temp = second->next;
        second->next = second->next->next;
        delete temp;
        return dum->next;
    }

    int numSquares(int n) {
        vector<int> dp(n + 1, INT_MAX);
        dp[0] = 0;
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j * j <= i; j++) {
                dp[i] = min(dp[i], dp[i - j * j] + 1);
            }
        }
        return dp[n];
    }

    int maximalRectangle(vector<vector<char>> &matrix) {

    }

    int search(vector<int> &nums, int target) {
        int l = 0, r = nums.size() - 1;
        while (l <= r) {
            int m = (l + r) / 2;
            if (nums[m] == target) {
                return m;
            }
            if (nums[l] <= nums[m]) {
                if (nums[l] <= target && target < nums[m]) {
                    r = m - 1;
                } else {
                    l = m + 1;
                }
            } else {
                if (nums[m] < target && target <= nums[r]) {
                    l = m + 1;
                } else {
                    r = m - 1;
                }
            }
        }
        return -1;
    }

    vector<int> twoSums(vector<int> &numbers, int target) {
        int n = numbers.size();
        int l = 0;
        int r = n - 1;
        while (l <= r) {
            if (numbers[l] + numbers[r] == target) {
                return {l + 1, r + 1};
            } else if (numbers[l] + numbers[r] > target) {
                r--;
            } else {
                l++;
            }
        }
        return {};
    }

    int removeDuplicates(vector<int> &nums) {
        int j = 1;
        for (int i = 1; i < nums.size(); i++) {
            if (nums[i] != nums[i - 1]) {
                nums[j] = nums[i];
                j++;
            }
        }
        return j;
    }

    int findMin(vector<int> &nums) {
        int n = nums.size();
        int l = 0, r = n - 1;
        while (l < r) {
            int m = l + (r - l) / 2;
            if (nums[m] >= nums[r]) {
                l = m + 1;
            } else {
                r = m;
            }
        }
        return nums[l];
    }

    int rangeSumBST(TreeNode *root, int low, int high) {
        if (!root) {
            return 0;
        }
        int right = rangeSumBST(root->right, low, high);
        int left = rangeSumBST(root->left, low, high);
        if (root->val >= low && root->val <= high) {
            return right + left + root->val;
        }
        return right + left;
    }

    vector<int> partitionLabels(string s) {
        map<char, int> last;
        vector<int> res;
        for (int i = 0; i < s.size(); i++) {
            last[s[i]] = i;
        }
        int edp = last[s[0]];
        int start = 0;
        while (start < s.size()) {
            edp = last[s[start]];
            int cnt = 0;
            while (start <= edp) {
                cnt++;
                edp = max(edp, last[s[start]]);
                start++;
            }
            res.push_back(cnt);
        }
        return res;
    }

    vector<int> productExceptSelf(vector<int> &nums) {
        int n = nums.size();
        vector<int> pref(n, 1);
        vector<int> suf(n, 1);
        for (int i = 1; i < n; i++) {
            pref[i] = pref[i - 1] * nums[i - 1];
        }
        for (int i = n - 2; i >= 0; i--) {
            suf[i] = suf[i + 1] * nums[i + 1];
        }
        vector<int> ans(n);
        for (int i = 0; i < n; i++) {
            ans[i] = pref[i] * suf[i];
        }
        return ans;
    }

    int maxPathSum(TreeNode *root) {
        int maxi = INT_MIN;
        path(maxi, root);
        return maxi;
    }

    int path(int &maxi, TreeNode *root) {
        if (!root) {
            return 0;
        }
        int l = max(path(maxi, root->left), 0);
        int r = max(path(maxi, root->right), 0);
        maxi = max(maxi, l + r + root->val);
        return root->val + max(l, r);
    }

    bool checkSubarraySum(vector<int> &nums, int k) {
        unordered_map<int, int> st;
        st[nums[0] % k] = 1;
        for (int i = 1; i < nums.size(); i++) {
            nums[i] += nums[i - 1];
            if (st[nums[i] % k] != 0) {
                if (abs(st[nums[i] % k] - (i + 1)) > 1) {
                    return 1;
                }

            } else {
                st[nums[i] % k] = i + 1;
            }
            if (nums[i] % k == 0) {
                return 1;
            }
        }
        return 0;
    }

    string reverseWords(string s) {
        string str = "";
        string check = "";
        int n = s.size();
        for (int i = 0; i < n; i++) {
            if (s[i] == ' ' || i == n - 1) {
                if (i == n - 1) {
                    check += s[i];
                }
                std::reverse(check.begin(), check.end());
                str += check;
                if (i != n - 1) {
                    str += ' ';
                }
                check = "";
            } else {
                check += s[i];
            }
        }
        return str;
    }

    string addStrings(string num1, string num2) {
        int n = num1.size() - 1;
        int m = num2.size() - 1;
        int carry = 0;
        string ans = "";
        while (n >= 0 || m >= 0) {
            int a = (n < 0) ? 0 : num1[n] - '0';
            int b = (m < 0) ? 0 : num2[m] - '0';
            int sum = a + b + carry;
            ans += to_string(sum % 10);
            carry = sum / 10;
            m--, n--;

        }
        if (carry > 0) {
            ans += to_string(carry);
        }
        std::reverse(ans.begin(), ans.end());
        return ans;
    }

    bool isPalindrome(ListNode *head) {
        ListNode* slow = head;
        ListNode* fast = head;
        while(fast->next != nullptr && fast->next->next != nullptr){
            slow = slow->next;
            fast = fast->next->next;
        }
        ListNode* newHead = reverse(slow->next);
        ListNode* first = head;
        ListNode* second = newHead;
        while(second != nullptr){
            if(first->val != second->val){
                return false;
            }
            first = first->next;
            second = second->next;
        }
        return true;
    }

    ListNode *reverse(ListNode *head) {
        if (head == nullptr || head->next == nullptr) {
            return head;
        }
        ListNode *newHead = reverse(head->next);
        ListNode *front = head->next;
        front->next = head;
        head->next = nullptr;
        return newHead;
    }

    int firstUniqChar(string s) {
        map<char, int> ss;
        for(auto x: s){
            ss[x]++;
        }
        for(int i = 0; i < s.size();i++){
            if(ss[s[i]]==1){
                return i;
            }
        }
        return -1;
    }

};

class NestedInteger {
public:
    // Return true if this NestedInteger holds a single integer, rather than a nested list.
    bool isInteger() const;

    // Return the single integer that this NestedInteger holds, if it holds a single integer
    // The result is undefined if this NestedInteger holds a nested list
    int getInteger() const;

    // Return the nested list that this NestedInteger holds, if it holds a nested list
    // The result is undefined if this NestedInteger holds a single integer
    const vector<NestedInteger> getList() const;
};


class NestedIterator {
    vector<int> flatList;
    int iterator;
public:
    void flattenList(vector<NestedInteger> nestedList) {
        for (NestedInteger ele: nestedList) {
            if (ele.isInteger()) {
                flatList.push_back(ele.getInteger());
            } else {
                flattenList(ele.getList());
            }
        }
    }

    NestedIterator(std::vector<NestedInteger> &nestedList) {
        iterator = 0;
        flattenList(nestedList);
    }

    int next() {
        if (hasNext()) {
            return flatList[iterator++];
        } else {
            return NULL;
        }
    }

    bool hasNext() {
        return iterator < flatList.size();
    }
};

class MyQueue {
private:
    stack<int> stack1, stack2;
public:
    MyQueue() {

    }

    void push(int x) {
        stack1.push(x);
    }

    void transferElements() {
        if (stack2.empty()) {
            while (!stack1.empty()) {
                stack2.push(stack1.top());
                stack1.pop();
            }
        }
    }

    int pop() {
        transferElements();
        int el = stack2.top();
        stack2.pop();
        return el;
    }

    int peek() {
        transferElements();
        return stack2.top();
    }

    bool empty() {
        return stack2.empty() && stack1.empty();
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
    vector<char> aaa = {'a', 'a', 'b', 'b', 'c', 'c', 'c'};
    Solution::compress(aaa);
}
