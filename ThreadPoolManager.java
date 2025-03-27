import java.util.concurrent.*;

public class ThreadPoolManager {
    private static int availableCores = Runtime.getRuntime().availableProcessors();
    private static int defaultWorkerThreads = availableCores; // 기본적으로 가용한 코어 수 사용

    private static ExecutorService executor;
    private static ExecutorService resultExecutor;

    // 스레드 풀 초기화 (자동 설정 or 사용자 지정)
    public static void initialize(int numThreads) {
        if (numThreads <= 0) numThreads = defaultWorkerThreads;

        executor = Executors.newFixedThreadPool(numThreads);

        // 설정된 스레드 개수 출력
        System.out.println("ThreadPoolManager | Worker Threads: " + numThreads);
    }

    // 일반 작업 스레드 풀
    public static ExecutorService getExecutor() {
        return executor;
    }

    // 결과 처리 스레드 풀
    public static ExecutorService getResultExecutor() {
        return resultExecutor;
    }

    // 스레드 풀 종료
    public static void shutdown() {
        executor.shutdown();
    }

    public static void awaitTermination() throws InterruptedException {
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }
}
