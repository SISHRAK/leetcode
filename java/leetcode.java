import java.util.*;

//  Definition for a binary tree node.
class TreeNode {
    int val;
    TreeNode left;
    TreeNode right;

    TreeNode() {
    }

    TreeNode(int val) {
        this.val = val;
    }

    TreeNode(int val, TreeNode left, TreeNode right) {
        this.val = val;
        this.left = left;
        this.right = right;
    }
}

class Solution {
    public boolean reorderedPowerOf2(int n) {
        char[] number = SortedNum(n);
        for (int i = 0; i < 30; i++) {
            char[] p = SortedNum(1 << i);
            if (Arrays.equals(number, p)) {
                return true;
            }
        }
        return false;
    }

    private char[] SortedNum(int n) {
        char[] digits = String.valueOf(n).toCharArray();
        Arrays.sort(digits);
        return digits;
    }

    public double[] convertTemperature(double celsius) {
        return new double[]{celsius + 273.15, celsius * 1.80 + 32.00};
    }

    public int numberOfCuts(int n) {
        if (n == 1) return 0;
        if (n % 2 == 0) return n / 2;
        return n;
    }

    public int[][] largestLocal(int[][] grid) {
        int n = grid.length;
        int res[][] = new int[n - 2][n - 2];
        for (int i = 0; i < n - 2; i++) {
            for (int j = 0; j < n - 2; j++) {
                res[i][j] = findMax(grid, i, j);
            }
        }
        return res;
    }

    public int findMax(int[][] grid, int i, int j) {
        int max = Integer.MIN_VALUE;
        for (int a = i; a < i + 3; a++) {
            for (int b = j; b < j + 3; b++) {
                max = Math.max(grid[a][b], max);
            }
        }
        return max;

    }

    public int commonFactors(int a, int b) {
        int cnt = 0;
        for (int i = 1; i <= Math.min(a, b); i++) {
            if (a % i == 0 && b % i == 0) {
                cnt++;
            }
        }
        return cnt;
    }

    public boolean evaluateTree(TreeNode root) {
        if (root.left == null && root.right == null) {
            return root.val == 0 ? false : true;
        } else if (root.val == 2) {
            return evaluateTree(root.left) || evaluateTree(root.right);
        }
        return evaluateTree(root.left) && evaluateTree(root.right);
    }

    public int mostFrequentEven(int[] nums) {
        int temp[] = new int[100001];
        for (int i : nums) {
            temp[i]++;
        }
        int min = 1000000, fr = 0, ans = -1;
        for (int i = 0; i < temp.length; i++) {
            if (i % 2 == 0 && fr < temp[i]) {
                fr = temp[i];
                ans = i;
            }
        }
        return ans;
    }

    public int maxDepth(TreeNode root) {
        if (root == null) return 0;
        int left = maxDepth(root.left);
        int right = maxDepth(root.right);
        return Math.max(left, right) + 1;
    }

    public int findTheDistanceValue(int[] arr1, int[] arr2, int d) {
        int ans = 0;
        for (int i = 0; i < arr1.length; i++) {
            for (int j = 0; j < arr2.length; j++) {
                if (Math.abs(arr1[i] - arr2[j]) <= d) {
                    ans++;
                }
            }
        }
        return arr1.length - ans;
    }

    public int countNumbersWithUniqueDigits(int n) {
        if (n == 0) {
            return 1;
        }
        if (n == 1) {
            return 10;
        }
        return (int) Math.pow(10, n) - 9;
    }

    public boolean increasingTriplet(int[] nums) {
        int a = Integer.MAX_VALUE;
        int b = Integer.MAX_VALUE;
        for (int i = 0; i < nums.length; i++) {
            if (nums[i] <= a) {
                a = nums[i];
            } else if (nums[i] < b) {
                b = nums[i];
            } else if (nums[i] > b) {
                return true;
            }
        }
        return false;
    }

    //    public int lengthOfLIS(int[] nums) {
//        int ans = 0;
//        int buf = 0;
//        for(int i = 0; i < nums.length;i++){
//            for(int j = i; j < nums.length;j++){
//                if()
//            }
//
//
//        }
//    }
    public int findDuplicate(int[] nums) {
        HashSet<Integer> set = new HashSet<>();
        for (int num : nums) {
            if (set.contains(num)) {
                return num;
            }
            set.add(num);
        }
        return -1;
    }

    public int integerBreak(int n) {
        if (n <= 1) {
            return 0;
        }
        if (n == 2) {
            return 1;
        }
        if (n == 3) {
            return 2;
        }
        int ans = 1;
        while (n > 4) {
            ans *= 3;
            n -= 3;
        }
        ans *= n;
        return ans;
    }

    public int minOperations(String[] logs) {
        List<String> f = new ArrayList<>();
        for (int i = 0; i < logs.length; i++) {
            if (logs[i].equals("./")) {
                continue;
            } else if (logs[i].equals("../")) {
                if (f.size() > 1) {
                    f.remove(f.size() - 1);
                } else {
                    f.clear();
                }
            } else {
                f.add(logs[i]);
            }
        }
        return f.size();
    }

    public int getMaximumGenerated(int n) {
        if (n == 0) return 0;
        if (n == 1) return 1;
        int[] a = new int[n + 1];
        a[0] = 0;
        a[1] = 1;
        int max = 0;
        for (int i = 2; i < n + 1; i++) {
            if (i % 2 == 0) {
                a[i] = a[i / 2];
            } else {
                a[i] = a[i / 2 + 1] + a[i / 2];
            }
            max = Math.max(max, a[i]);
        }
        return max;
    }

    public boolean arrayStringsAreEqual(String[] word1, String[] word2) {
        String a = "";
        String b = "";
        for (int i = 0; i < word1.length; i++) {
            a += word1[i];
        }
        for (int i = 0; i < word2.length; i++) {
            b += word2[i];
        }
        if (Arrays.equals(a.toCharArray(), b.toCharArray())) {
            return true;
        }
        return false;
    }

    public int maximumScore(int[] nums, int k) {
        int res = nums[k], mini = nums[k], i = k, j = k, n = nums.length;
        while (i > 0 || j < n - 1) {
            if (i == 0) {
                j++;
            } else if (j == n - 1) {
                i--;
            } else if (nums[i - 1] < nums[j + 1]) {
                j++;
            } else {
                i--;
            }
            mini = Math.min(mini, Math.min(nums[i], nums[j]));
            res = Math.max(res, mini * (j - i + 1));
        }
        return res;
    }

    public int xorBeauty(int[] nums) {
        int ans = 0;
        for (int num : nums) {
            ans ^= num;
        }
        return ans;
    }

    public int minDistance(int[] houses, int k) {
        Arrays.sort(houses);
        int n = houses.length;
        int[] dp = new int[n];
        for (int i = 1; i < n; i++) {
            dp[i] = dp[i - 1] + houses[i] - houses[i / 2];
        }
        for (int i = 0; i < k - 1; i++) {
            int[] next = new int[n];
            for (int j = 0; j < next.length; j++) {
                next[j] = Integer.MAX_VALUE;
            }
            for (int j = 0; j < n; j++) {
                int sum = 0;
                for (int m = j; m >= 0; m--) {
                    sum += houses[(m + j + 1) >> 1] - houses[m];
                    int buf;
                    if (m == 0) {
                        buf = 0;
                    } else {
                        buf = dp[m - 1];
                    }
                    next[j] = Math.min(next[j], buf + sum);
                }
            }
            dp = next;
        }
        return dp[n - 1];
    }

    public int[] shortestToChar(String s, char c) {
        int[] ans = new int[s.length()];
        char[] aa = s.toCharArray();
        int buf = 0;
        int prev = s.length();
        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == c) {
                prev = 0;
                ans[i] = 0;
            } else {
                ans[i] = ++prev;
            }
        }
        prev = s.length();
        for (int i = s.length() - 1; i >= 0; i--) {
            if (s.charAt(i) == c) {
                prev = 0;
                ans[i] = Math.min(ans[i], 0);
            } else
                ans[i] = Math.min(ans[i], ++prev);
        }
        return ans;
    }


    public int numBusesToDestination(int[][] routes, int source, int target) {
        Map<Integer, List<Integer>> stopToBuses = new HashMap<>();
        for (int busId = 0; busId < routes.length; busId++) {
            for (int stop : routes[busId]) {
                stopToBuses.computeIfAbsent(stop, k -> new ArrayList<>()).add(busId);
            }
        }
        if (!stopToBuses.containsKey(source) || !stopToBuses.containsKey(target)) {
            return -1;
        }
        if (source == target) {
            return 0;
        }
        Queue<Integer> queue = new LinkedList<>();
        Set<Integer> busesTaken = new HashSet<>();
        Set<Integer> stopsVisited = new HashSet<>();
        int res = 0;
        queue.offer(source);
        while (!queue.isEmpty()) {
            res++;
            int stopsToProcess = queue.size();
            for (int i = 0; i < stopsToProcess; i++) {
                int currentStop = queue.poll();
                for (int busId : stopToBuses.getOrDefault(currentStop, new ArrayList<>())) {
                    if (busesTaken.contains(busId)) {
                        continue;
                    }
                    busesTaken.add(busId);
                    for (int nextStop : routes[busId]) {
                        if (stopsVisited.contains(nextStop)) {
                            continue;
                        }
                        if (nextStop == target) {
                            return res;
                        }
                        queue.offer(nextStop);
                        stopsVisited.add(nextStop);
                    }
                }
            }
        }
        return -1;
    }

    public int numMatchingSubseq(String s, String[] words) {
        Map<String, Integer> map = new HashMap<>();
        int ans = 0;
        for (String word : words) {
            map.put(word, map.getOrDefault(word, 0) + 1);
        }
        for (String word : map.keySet()) {
            int i = 0, j = 0;
            while (i < s.length() && j < word.length()) {
                if (s.charAt(i) == word.charAt(j)) {
                    j++;
                }
                i++;
            }
            if (j == word.length()) {
                ans += map.get(word);
            }
        }
        return ans;
    }
}
