class Solution {
    fun stringMatching(words: Array<String>): List<String> {
        return words.sorted().filter{word -> words.count{ it.contains(word) } > 1}
    }
    fun canBeEqual(target: IntArray, arr: IntArray): Boolean{
        return target.sorted() == arr.sorted();
    }
    fun maxSatisfaction(satisfaction: IntArray): Int {
        val sorted = satisfaction.sortedDescending()
        var level = 0
        var sum = 0
        var res = 0
        for(x in sorted){
            level += x
            sum += level
            res = maxOf(sum, res)
            if(level < 0){
                break
            }
        }
        return res
    }
}

fun main(args: Array<String>) {

}
