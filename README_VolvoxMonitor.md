# 🌿 Volvox Monitor — Система контроля параметров биореактора

> **Платформа:** Arduino UNO (тестирование) → NodeMCU LoLin V3 ESP8266 (постоянная работа + WiFi)
> **Датчики:** DS18B20 (температура) · MH-Z19B (CO₂) · pH модуль · Фоторезистор
> **Данные:** ThingSpeak → онлайн-графики в браузере

---

## 📦 Необходимые библиотеки

Установить через **Arduino IDE → Инструменты → Управлять библиотеками**:

| Библиотека | Автор | Для чего |
|---|---|---|
| `OneWire` | Jim Studt | Протокол 1-Wire для DS18B20 |
| `DallasTemperature` | Miles Burton | Чтение DS18B20 |
| `MHZ19` | Jonathan Dempsey | Чтение MH-Z19B |
| `ThingSpeak` | MathWorks | Отправка данных в облако |
| `SoftwareSerial` | Встроена | Второй UART для CO₂ |

Для NodeMCU дополнительно:
**Arduino IDE → Файл → Настройки → Доп. ссылки для менеджера плат:**
```
https://arduino.esp8266.com/stable/package_esp8266com_index.json
```
Затем: **Инструменты → Плата → Менеджер плат** → найти `esp8266` → Установить

---

# ШАГ 1 — Тестирование на Arduino UNO

> Тестируем каждый датчик **по отдельности** на Arduino UNO.
> Результаты видны в **Инструменты → Монитор порта → 9600 baud**.

---

## ТЕСТ 1.1 — DS18B20 (температура)

### Подключение

```
DS18B20             Arduino UNO
──────────────────────────────────
Красный  (VCC)  →   5V
Чёрный   (GND)  →   GND
Жёлтый   (DATA) →   D2

Резистор 4.7 кОм: между 5V и D2 (pull-up, обязательно!)
```

```
  5V ──┬──[ 4.7kΩ ]──┬── D2 (UNO)
       │             │
      VCC           DATA     DS18B20
       │             │
      GND ──────── GND
```

### Код

```cpp
#include <OneWire.h>
#include <DallasTemperature.h>

#define SENSOR_PIN 2          // Жёлтый провод → D2

OneWire oneWire(SENSOR_PIN);
DallasTemperature ds(&oneWire);

void setup() {
  Serial.begin(9600);
  ds.begin();
  Serial.println("DS18B20 готов");
}

void loop() {
  ds.requestTemperatures();
  float t = ds.getTempCByIndex(0);
  Serial.print("Температура: ");
  Serial.print(t);
  Serial.println(" °C");
  delay(1000);
}
```

### Ожидаемый результат
```
DS18B20 готов
Температура: 22.50 °C
Температура: 22.56 °C
```

---

## ТЕСТ 1.2 — MH-Z19B (CO₂)

### Подключение

```
MH-Z19B             Arduino UNO
──────────────────────────────────
VCC (5V)  →   5V
GND       →   GND
TX        →   D2   (SoftwareSerial RX)
RX        →   D3   (SoftwareSerial TX)

⚠️ Пины D0/D1 НЕ использовать — они заняты USB
```

### Код

```cpp
#include <SoftwareSerial.h>
#include <MHZ19.h>

SoftwareSerial co2Serial(2, 3);  // RX=D2, TX=D3
MHZ19 mhz;

void setup() {
  Serial.begin(9600);
  co2Serial.begin(9600);
  mhz.begin(co2Serial);
  mhz.autoCalibration(false);   // ⚠️ выкл. для биореактора!
  Serial.println("CO2 готов. Прогрев 30 сек...");
  delay(30000);
}

void loop() {
  int co2 = mhz.getCO2();
  Serial.print("CO2: ");
  Serial.print(co2);
  Serial.println(" ppm");
  delay(2000);
}
```

### Ожидаемый результат
```
CO2 готов. Прогрев 30 сек...
CO2: 412 ppm
CO2: 415 ppm
```
> 400–450 ppm = норма воздуха. После помещения в биореактор значения изменятся.

---

## ТЕСТ 1.3 — pH датчик (аналоговый)

### Подготовка электрода (ВАЖНО, сделать за день до теста!)
```
1. Погрузить pH-электрод в 1% раствор HCl на 24 часа → активация
2. Промыть дистиллированной водой
3. Откалибровать в буферных растворах pH 4.00 и pH 7.00
```

### Подключение

```
pH модуль           Arduino UNO
──────────────────────────────────
VCC       →   5V
GND       →   GND
Signal    →   A0

⚠️ На UNO делитель НЕ нужен — A0 читает 0–5V напрямую
⚠️ На NodeMCU делитель НУЖЕН (см. Шаг 2)
```

### Код

```cpp
void setup() {
  Serial.begin(9600);
}

void loop() {
  int raw = analogRead(A0);                    // 0–1023
  float voltage = raw * (5.0 / 1023.0);        // → в вольты

  // Калибровочная формула (уточнить после калибровки):
  // pH 7  ≈ 2.50 V
  // pH 4  ≈ 3.05 V  → наклон ≈ 0.18 V/pH
  float pH = 7.0 + ((2.5 - voltage) / 0.18);

  Serial.print("Сырое: "); Serial.print(raw);
  Serial.print("  Напряжение: "); Serial.print(voltage, 3); Serial.print(" V");
  Serial.print("  pH: "); Serial.println(pH, 2);
  delay(1000);
}
```

### Ожидаемый результат
```
Сырое: 512  Напряжение: 2.501 V  pH: 7.00
Сырое: 491  Напряжение: 2.398 V  pH: 7.57
```

---

# ШАГ 2 — NodeMCU + WiFi + онлайн-графики

## Настройка ThingSpeak

```
1. Перейти на https://thingspeak.mathworks.com
2. Зарегистрироваться (бесплатно)
3. New Channel → заполнить:
   Field 1: Temperature (°C)
   Field 2: CO2 (ppm)
   Field 3: pH
   Field 4: Light
4. Сохранить → API Keys → скопировать Write API Key
5. Запомнить Channel ID (число в URL)
```

---

## Подключение всех датчиков к NodeMCU

```
DS18B20 (температура)
──────────────────────────────────────────
Красный  (VCC)  →  NodeMCU  3.3V
Чёрный   (GND)  →  NodeMCU  GND
Жёлтый   (DATA) →  NodeMCU  D1  (GPIO5)
Резистор 4.7кОм: между 3.3V и D1

MH-Z19B (CO₂)
──────────────────────────────────────────
VCC (5V)  →  NodeMCU  VIN   ← 5V из USB!
GND       →  NodeMCU  GND
TX        →  NodeMCU  D6   (GPIO12)  RX-вход NodeMCU
RX        →  NodeMCU  D5   (GPIO14)  TX-выход NodeMCU
✅ LLC-преобразователь НЕ нужен — MH-Z19B совместим с 3.3V логикой

pH датчик (аналоговый) + ДЕЛИТЕЛЬ НАПРЯЖЕНИЯ
──────────────────────────────────────────
VCC       →  NodeMCU  VIN  (5V)
GND       →  NodeMCU  GND
Signal    →  R1 (1кОм) → A0 → R2 (2кОм) → GND

Схема делителя:
  Signal ──[ R1=1kΩ ]──┬── A0 (NodeMCU)
                       │
                    [ R2=2kΩ ]
                       │
                      GND

⚠️ Зачем делитель: pH даёт до 3.5V, A0 NodeMCU читает макс. 3.3V
   Без делителя: показания врут, возможно повреждение ADC
   С делителем: 3.5V × 2/(1+2) = 2.33V → безопасно ✅
   R1=1кОм и R2=2×1кОм есть в наборе Arduino UNO!

Фоторезистор 5516
──────────────────────────────────────────
Нога 1    →  NodeMCU  3.3V
Нога 2    →  NodeMCU  D7  (GPIO13)
10кОм     →  между D7 и GND

Питание NodeMCU
──────────────────────────────────────────
microUSB  →  Зарядка 5V / 2A (от розетки 220V)
⚠️ Минимум 1A, рекомендуется 2A — датчики добавляют нагрузку
```

---

## Полный код для NodeMCU

```cpp
// ╔══════════════════════════════════════════════════════════╗
// ║         VOLVOX MONITOR — NodeMCU + ThingSpeak           ║
// ║  Датчики: DS18B20 · MH-Z19B · pH аналог · Фоторезистор ║
// ╚══════════════════════════════════════════════════════════╝

#include <ESP8266WiFi.h>
#include <ThingSpeak.h>
#include <OneWire.h>
#include <DallasTemperature.h>
#include <SoftwareSerial.h>
#include <MHZ19.h>

// ── НАСТРОЙКИ WiFi и ThingSpeak ────────────────────────────
const char* WIFI_SSID  = "ИМЯ_ВАШЕЙ_СЕТИ";   // ← заменить
const char* WIFI_PASS  = "ПАРОЛЬ_WIFI";       // ← заменить
unsigned long CHANNEL_ID = 123456;            // ← ваш Channel ID
const char*   API_KEY    = "XXXXXXXXXXXX";    // ← Write API Key

// ── ПИНЫ NodeMCU ───────────────────────────────────────────
#define DS_PIN    D1   // DS18B20: жёлтый провод + 4.7кОм к 3.3V
#define CO2_RX    D6   // MH-Z19B TX → D6 (SoftwareSerial RX)
#define CO2_TX    D5   // MH-Z19B RX ← D5 (SoftwareSerial TX)
#define PHOTO_PIN D7   // Фоторезистор + 10кОм к GND
// A0 ← pH Signal (через делитель R1=1кОм, R2=2кОм)

// ── ОБЪЕКТЫ ────────────────────────────────────────────────
OneWire          oneWire(DS_PIN);
DallasTemperature ds(&oneWire);
SoftwareSerial   co2Ser(CO2_RX, CO2_TX);
MHZ19            mhz;
WiFiClient       client;

// ── SETUP ──────────────────────────────────────────────────
void setup() {
  Serial.begin(9600);
  Serial.println("
=== Volvox Monitor ===");

  // Инициализация датчиков
  ds.begin();
  co2Ser.begin(9600);
  mhz.begin(co2Ser);
  mhz.autoCalibration(false);   // ⚠️ отключить для замкнутого биореактора!
  pinMode(PHOTO_PIN, INPUT);

  // Подключение к WiFi
  WiFi.begin(WIFI_SSID, WIFI_PASS);
  Serial.print("Подключение к WiFi");
  while (WiFi.status() != WL_CONNECTED) {
    delay(500);
    Serial.print(".");
  }
  Serial.println();
  Serial.print("✅ WiFi подключён. IP: ");
  Serial.println(WiFi.localIP());

  // Инициализация ThingSpeak
  ThingSpeak.begin(client);

  // Прогрев CO₂ датчика
  Serial.println("⏳ Прогрев MH-Z19B (30 сек)...");
  delay(30000);
  Serial.println("✅ Готов к работе!");
}

// ── LOOP ───────────────────────────────────────────────────
void loop() {

  // 1. Температура (DS18B20)
  ds.requestTemperatures();
  float temp = ds.getTempCByIndex(0);

  // 2. CO₂ (MH-Z19B)
  int co2 = mhz.getCO2();

  // 3. pH (аналоговый, через делитель)
  int   phRaw   = analogRead(A0);
  // NodeMCU A0: 0–1023 → 0–3.3V (LoLin v3 имеет встроенный делитель)
  // Компенсируем внешний делитель R1/(R1+R2) = 1/3 → умножаем на 1.5
  float voltage = phRaw * (3.3 / 1023.0) * 1.5;
  float pH      = 7.0 + ((2.5 - voltage) / 0.18);  // уточнить после калибровки

  // 4. Освещённость (фоторезистор, цифровой порог)
  int light = digitalRead(PHOTO_PIN);  // 1=светло, 0=темно

  // Вывод в Serial Monitor
  Serial.printf("🌡 T=%.1f°C  🌬 CO2=%d ppm  🧪 pH=%.2f  💡 Свет=%s
",
                temp, co2, pH, light ? "ДА" : "НЕТ");

  // Отправка в ThingSpeak
  ThingSpeak.setField(1, temp);
  ThingSpeak.setField(2, co2);
  ThingSpeak.setField(3, pH);
  ThingSpeak.setField(4, light);
  int httpCode = ThingSpeak.writeFields(CHANNEL_ID, API_KEY);

  if (httpCode == 200) {
    Serial.println("   → ✅ Данные отправлены в ThingSpeak");
  } else {
    Serial.print("   → ❌ Ошибка отправки. HTTP код: ");
    Serial.println(httpCode);
  }

  delay(30000);  // ⚠️ ThingSpeak: минимальный интервал 15 сек
}
```

---

## Как смотреть графики

```
1. thingspeak.com → My Channels → ваш канал → Charts
   Автоматические графики по каждому полю ✅

2. Телефон: приложение ThingView (Android / iOS, бесплатно)
   Показывает все 4 графика в реальном времени

3. Экспорт данных: Export → CSV → открыть в Excel/Python
```

---

## Устранение типичных проблем

| Проблема | Причина | Решение |
|---|---|---|
| `DS18B20 показывает -127°C` | Нет резистора 4.7кОм | Проверить pull-up между VCC и DATA |
| `CO2 = 0 или -1` | Прогрев не завершён | Подождать 60–90 сек после включения |
| `pH показывает мусор` | Электрод не активирован | 24 часа в 1% HCl, затем калибровка |
| `WiFi не подключается` | Неверный пароль | Проверить WIFI_SSID и WIFI_PASS |
| `HTTP 0 или -301` | Нет интернета | Проверить WiFi на другом устройстве |
| `ThingSpeak ошибка 0` | Слишком частая отправка | delay() должен быть ≥ 15000 мс |

---

*Volvox Monitor v1.0 | Сколтех | 2026*
