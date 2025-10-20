# Student information

---
 - **Suchan Roh** 
 - 2022094839 `id`
 - 졸업프로젝트2025-2팀
 - 졸업프로젝트 지도 교수: **백은옥 교수님**

# Project information

---
## Graduation project


- Duration: `2025-1` ~ `2025-2`
- Project link: 
- https://github.com/tncks/Capstone_MODPlus
- Owner (License): 
  - This graduation project forked from:
  - https://github.com/HanyangBISLab/MODplus
  - **나승진 박사님**
  
---

# Full documentation

## 1. 개요

### 기본 정보
- **프로젝트명**: MODPlus
- **생성일**: 2025년 3월 24일
- **최종 업데이트**: 2025년 10월 20일


### 기술 스택
- **언어**: Java
- **Java 버전**: Java 23
- **빌드 도구**: Java compiler, `intelliJ`
- **외부 라이브러리**:
    - `jdom-1.1.3.jar` 
    - `jrap_StAX_v5.2.jar` 
- **동시성 프레임워크**: `java.util.concurrent` (`ThreadPoolExecutor`)


---


## 2. 프로젝트 수행 결과

### 주요 개선점
- **멀티스레드 아키텍처**: CPU 코어 활용 극대화  
- **정확도 100% 기반 일관성 유지**: 결과 출력 순서 보장   
- **스레드 안전성**: ThreadLocal 패턴 적용

### 동시성 처리
- **동시성 프로그래밍**: ConcurrentHashMap, AtomicInteger 활용
- **메모리 최적화**: ThreadLocal 정리로 유휴 자원 반납
- **I/O 최적화**: 비동기 파일 쓰기 (종료 직전에 모아 일괄 쓰기 처리)



---

## 3. 기존 대비 업데이트 사항

### 주요 변경 내역

#### **초기 상태 (기존)**
- 단일 스레드 사용
- 전역 `static` 변수로 상태 관리
- 순차적 스펙트럼 분석

#### **중간 상태**
- **멀티스레드 도입**: `ThreadPoolExecutor` 기반
- **TaskDecomposition**: `ModPlusTask` 클래스로 작업 분리
- **결과 수집**: `ConcurrentHashMap` 사용
- **I/O 최적화**: CPU 작업 후 일괄 쓰기 처리

#### **현재 상태 (Today)**
- **스레드 안전성**: `ThreadLocal` 패턴
- **Null 안전 처리**: 예외 처리 패턴
- **정확도 100%**: 검증된 알고리즘
- **모듈화**: 함수 모듈화

---

## 4. 메인 함수의 작동 방식

**주요 변경사항 (`main` 리팩토링)**:
```java
// (MODPlus.java) 단일 스레드 → 멀티스레드 전환
ThreadPoolExecutor executor = new ThreadPoolExecutor(
    corePoolSize, maxPoolSize, keepAliveTime, TimeUnit.SECONDS,
    new LinkedBlockingQueue<>(qCapacity)
);

while (scaniter.hasNext()) {
    ArrayList<MSMScan> chargedSpectra = scaniter.getNext();
    executor.submit(new ModPlusTask(chargedSpectra, ixPDB, considerIsotopeErr, scaniter.getIndex()));
}
```

---

**기타 변경 사항**:
- **ThreadLocal 패턴 도입**: `Mutables` 클래스를 스레드 로컬 변수로 관리하여 경합 조건 해결
- **ConcurrentHashMap 사용**: 스레드 안전한 결과 저장
- **성능 향상**: 병렬 처리 적용
- **I/O 분리**: 별도 스레드로 파일 쓰기 작업 처리 (마지막 단계 도달 시 일괄 처리)

**변경 파일**:
- `MODPlus.java` 
- 새로운 클래스 추가: `ModPlusTask`, `ResultEntry`

---

#### **정확도 100% 달성 (2025-05-21) - PTMDB 알고리즘 개선**


**개선사항**:
```java
// (PTMDB.java) PTM 검색 로직 리팩토링
public PTMSearchResult searchPTM(Sequence seq, double massDiff, PTMPosition position) {
    // Hash 캐싱 제거 → 스레드 안전성 강화
    PTMSearchResult searchResult;
    ArrayList<PTMRun> newGapInterpret = new ArrayList<>();
    
    double error = Math.abs(massDiff);
    if (error < (ThreadLocalMutables.get().nonModifiedDelta)) {
        return new PTMSearchResult(newGapInterpret, true);
    }
    
    // DFS 탐색 후 결과 즉시 반환
    context.findPTM_DFS(0.0, 0, 0, 0);
    context.setResetResultToNull(); // 메모리 누수 방지
    
    threadLocalContext.remove(); // ThreadLocal 정리
    return searchResult;
}
```

**추가된 파일**:
- `GapType.java` (Enum 타입 정의)
- `DeepCloneObj.java` (객체 복제 유틸리티)
- `ProteinContainer.java` (단백질 관련 항목)


**알고리즘 변경**:
- **MGFIterator 처리**: 일부 변경 (스레딩 기반)
- **PTMDB 리팩토링**: 일부 변경 (캐싱 전략 변경)

---

#### **기타 개선 사항**

| 커밋 | 날짜 | 내용 |
|------|------|------|
| `ec720d1` | 2025-05-21 | `run()` 함수 모듈화 및 함수 추가 |
| `82b33da` | 2025-04-13 | Print 문 개선 (로깅 가독성) |
| `a54d5ff` | 2025-04-13 | 작은 버그 수정 (chore) |
| `3ffbb31` | 2025-04-13 | 정확도 이슈 해결 (Accuracy issue resolved) |
| `39750525` | 2025-04-01 | 설명 텍스트 추가 (문서화) |

---

## 5. 트러블슈팅 (각종 버그 해결)

### 실행 시 산출된 결과물의 정확도 저하 이슈

| 발생 시기 | 원인        | 해결 방법 |
|-----------|-----------|-----------|
| 2025-04-13 | 스펙트럼 파싱 오류 | MGFIterator 로직 수정 |
| 2025-05-21 | PTM 검색 캐싱 문제 | 캐싱 전략 제거 + ThreadLocal |

**세부 원인 식별 경로**:
- **멀티스레드 환경**에서 `static` 변수 공유로 인한 경합 조건 발견 (스레드 덤프 툴 사용)
- **ThreadLocal** 도입 전 스레드 간 상태 오염 (소수의 candidates 가 누락되었던 문제)

**최종 해결책**:
```java
// (MODPlus.java) ThreadLocalMutables 패턴
Mutables mutables = ThreadLocalMutables.get();
mutables.maxNoOfC13 = chargedSpectra.get(z).getMaxNoOfC13();
// ... 스레드별 독립적 상태 관리

// finally 블록에서 반드시 정리
finally {
    ThreadLocalMutables.clear(); // 메모리 누수 방지
}
```

**이슈 해결 일자**:
- 2025-08-13  `06a6600`
---

### 쓰기 작업 성능 최적화

**I/O 처리 개선**:
- "timeout" (타임아웃 처리 개선)
- "I/O handle" (I/O 작업 분리)

**개선 사항**:
- **CPU 작업과 I/O 작업 분리**:
   ```java
   ScheduledExecutorService scheduler = Executors.newSingleThreadScheduledExecutor();
   scheduler.schedule(ioHandler, 1_000, TimeUnit.MILLISECONDS); // 1초 지연 후 I/O
   ```



---



## 6. 성능 측정



**개선 전 (단일 스레드)**:
```java
// (MODPlus.java) 순차 처리
while (scaniter.hasNext()) {
    // 각 스캔을 순차적으로 처리
    
}
```

**개선 후 (멀티 스레드)**:
```java
// (MODPlus.java) 병렬 처리
ThreadPoolExecutor executor = new ThreadPoolExecutor(
    Runtime.getRuntime().availableProcessors() - 1, ...
);

while (scaniter.hasNext()) {
    executor.submit(new ModPlusTask(...));
}
// 시간: CPU 코어 개수에 근접해 비례하는 수준의 실행 시간 개선
```

**성능 향상률**: 최대 `27.0 배` (단, `36` 코어 환경)

**평균 성능 향상률**: 일반적으로 약 `25.0배 ~ 26.0배` 사이에 분포


---

## 7. 프로젝트 진행 절차 

### a. 프로파일링 진행
### 리팩토링을 진행하기 위해, 코드 병목 지점을 먼저 식별.
#### `intelliJ Runtime Profiler`
#### 프로파일러를 통해 평균 실행 시간보다 느린 속도를 보이는 함수 식별에 성공.
- #### 위 분석을 바탕으로 `MODPlus.java` 파일의 병목 지점을 리팩토링하기로 계획.

### b. 구현
### 식별된 병목 지점에 대해 실제 리팩토링을 진행.
#### 예:
```java
// (MODPlus.java) 병렬 처리
ThreadPoolExecutor executor = new ThreadPoolExecutor(
    Runtime.getRuntime().availableProcessors() - 1, ...
);

while (scaniter.hasNext()) {
    executor.submit(new ModPlusTask(...));
}
// 예상 시간: N개 스캔 × 평균 처리 시간 / (코어 수 - 1)
```



---

## 8. Appendix

### 외부 라이브러리 (변경사항 없이, 그대로 사용 유지)

| 라이브러리 | 버전 | 용도 | 크기 |
|-----------|------|------|------|
| JDOM | 1.1.3 | XML 파싱 (param.xml) | 151 KB |
| jRAP StAX | 5.2 | MZXML 파일 파싱 | 157 KB |

### Java 표준 라이브러리 (변경점)

**사용된 동시성 패키지**:
```java
import java.util.concurrent.*;           // ThreadPoolExecutor, ConcurrentHashMap
import java.util.concurrent.atomic.*;     // AtomicInteger
```

**동시성 패턴**:
- `ThreadPoolExecutor` → 작업 큐 관리
- `ConcurrentHashMap` → 스레드 기반 일관된 결과 저장
- `AtomicInteger` → 스레드 안전 카운터
- `ThreadLocal` → 스레드별 상태 분리
