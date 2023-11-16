import java.util.Arrays;
import java.util.HashSet;

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
        for(int i = 0; i < arr1.length;i++){
            for(int j = 0; j < arr2.length;j++){
                if(Math.abs(arr1[i] - arr2[j]) <= d){
                    ans++;
                }
            }
        }
        return arr1.length - ans;
    }
    public int countNumbersWithUniqueDigits(int n) {
        if (n == 0){
            return 1;
        }
        if(n == 1){
            return 10;
        }
        return (int) Math.pow(10, n) - 9;
    }
    public boolean increasingTriplet(int[] nums) {
        int a = Integer.MAX_VALUE;
        int b = Integer.MAX_VALUE;
        for(int i = 0; i < nums.length;i++){
            if(nums[i] <= a){
                a = nums[i];
            }
            else if(nums[i] < b){
                b = nums[i];
            }
            else if(nums[i] > b){
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
        for(int num : nums){
            if(set.contains(num)){
                return num;
            }
            set.add(num);
        }
        return -1;
    }
    public int integerBreak(int n) {
        if(n <= 1){
            return 0;
        }
        if (n == 2) {
            return 1;
        }
        if(n == 3){
            return 2;
        }
        int ans = 1;
        while(n > 4){
            ans *= 3;
            n -= 3;
        }
        ans *= n;
        return ans;
    }
}
