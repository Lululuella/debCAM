import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
//import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
//import java.util.Random;
 
public class FixSizedPriorityQueue<E extends Comparable<? super E>> {
    private PriorityQueue<E> queue;
    private int maxSize; 
 
    public FixSizedPriorityQueue(int maxSize) {
        if (maxSize <= 0)
            throw new IllegalArgumentException();
        this.maxSize = maxSize;
        this.queue = new PriorityQueue<E>(maxSize, new Comparator<E>() {
            public int compare(E o1, E o2) {
                // max heap: o2-o1
                return (o2.compareTo(o1));
            }
        });
    }
 
    public void add(E e) {
        if (queue.size() < maxSize) { // not full, add directly
            queue.add(e);
        } else { // full
            E peek = queue.peek();
            if (e.compareTo(peek) < 0) { // max heap, keep smaller one
                queue.poll();
                queue.add(e);
            }
        }
    }
 
    public List<E> sortedList() {
        List<E> list = new ArrayList<E>(queue);
        Collections.sort(list); // PriorityQueue itself is not sorted
        return list;
    }
 
    /*
    public static void main(String[] args) {
        final FixSizedPriorityQueue<Integer> pq = new FixSizedPriorityQueue<Integer>(10);
        Random random = new Random();
        int rNum = 0;
        System.out.println("100 rnd:-----------------------------------");
        for (int i = 1; i <= 100; i++) {
            rNum = random.nextInt(1000);
            System.out.println(rNum);
            pq.add(rNum);
        }
        System.out.println("PriorityQueue:-----------------------------------");
        Iterable<Integer> iter = new Iterable<Integer>() {
            public Iterator<Integer> iterator() {
                return pq.queue.iterator();
            }
        };
        for (Integer item : iter) {
            System.out.print(item + ", ");
        }
        System.out.println();
        System.out.println("PriorityQueue sorted:-----------------------------------");
        
        //for (Integer item : pq.sortedList()) { System.out.println(item); }
        
        while (!pq.queue.isEmpty()) {
            System.out.print(pq.queue.poll() + ", ");
        }
    }
    */
}
